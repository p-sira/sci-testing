"""Fast, mock-compiler tests for interprocess native build caching.

All compilers and CDLL loads in this suite are mocked: never compile the
repository's expensive 40+ native translation units.
"""
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest import mock

from sci_testing import native_compare as native


def _fake_elf_run(command, *, events=None, delay=0., fail=False, **kwargs):
    """Stand-in for one invocation of the compiler, not a real ELF linker."""
    output = Path(command[command.index('-o') + 1])
    if events is not None:
        with open(events, 'a', encoding='ascii') as stream:
            stream.write('compile\n')
    if delay:
        time.sleep(delay)
    if fail:
        raise subprocess.CalledProcessError(1, command)
    output.write_bytes(b'\x7fELFmocked-native-library')
    return subprocess.CompletedProcess(command, 0)


def _concurrent_caller(root, queue, events):
    # On Linux fork and spawn both get an independently opened flock handle.
    try:
        with (mock.patch.object(native, 'ROOT', Path(root)),
              mock.patch.object(native, '_NATIVE_SOURCES', ('core.cpp',)),
              mock.patch.object(native, '_compiler_command',
                                return_value=('/fake/c++',)),
              mock.patch.object(native, '_compiler_and_libraries',
                                return_value={'boost': '1_83', 'gsl': '2.7'}),
              mock.patch.object(native.subprocess, 'run',
                                side_effect=lambda args, **kwargs:
                                _fake_elf_run(args, events=events, delay=.25,
                                              **kwargs)),
              mock.patch.object(native.ct, 'CDLL',
                                side_effect=lambda path: str(path))):
            queue.put(('ok', native._native_library()))
    except BaseException as error:
        queue.put(('error', repr(error)))


class NativeBuildCacheTests(unittest.TestCase):
    def setUp(self):
        workspace = tempfile.TemporaryDirectory()
        self.addCleanup(workspace.cleanup)
        self.root = Path(workspace.name)
        (self.root / 'native').mkdir()
        self.source = self.root / 'native/core.cpp'
        self.source.write_text('int version() { return 1; }\n')
        self.run = mock.Mock(side_effect=_fake_elf_run)
        self.probe = mock.Mock(return_value={
            'compiler_version': 'g++ mock',
            'library_headers': 'Boost 1_83, GSL 2.7',
            'linked_libraries': [],
        })
        for item in (
            mock.patch.object(native, 'ROOT', self.root),
            mock.patch.object(native, '_NATIVE_SOURCES', ('core.cpp',)),
            mock.patch.object(native, '_compiler_command',
                              return_value=('/fake/c++',)),
            mock.patch.object(native, '_compiler_and_libraries', self.probe),
            mock.patch.object(native.subprocess, 'run', self.run),
            mock.patch.object(native.ct, 'CDLL',
                              side_effect=lambda path: str(path)),
        ):
            item.start()
            self.addCleanup(item.stop)

    def test_cache_hit_and_fingerprint_specific_binary(self):
        first = native._native_library()
        self.assertEqual(self.run.call_count, 1)
        self.assertTrue(Path(first).exists())
        self.assertTrue(first.endswith('.so'))
        self.assertIn('native-comparisons-', first)
        self.assertEqual(native._native_library(), first)
        self.assertEqual(self.run.call_count, 1)
        manifest = Path(first).with_suffix('.json')
        meta = json.loads(manifest.read_text())
        self.assertEqual(meta['sha256'], native._file_sha256(Path(first)))
        self.assertEqual(len(meta['fingerprint']), 64)
        arguments = self.run.call_args.args[0]
        self.assertEqual(arguments[:5],
                         ['/fake/c++', *native._NATIVE_FLAGS])
        self.assertEqual(arguments[-4:-2], list(native._NATIVE_LIBS))

        self.source.write_text('int version() { return 2; }\n')
        second = native._native_library()
        self.assertNotEqual(first, second)
        self.assertEqual(self.run.call_count, 2)
        self.assertTrue(Path(first).is_file(), 'loaded old inode remains intact')
        self.assertEqual(native._native_library(), second)
        self.assertEqual(self.run.call_count, 2)

    def test_version_compiler_and_flags_invalidate_cache(self):
        original = native._native_library()
        self.probe.return_value = {
            'compiler_version': 'g++ mock',
            'library_headers': 'Boost 1_83, GSL 2.8',
            'linked_libraries': [],
        }
        changed_headers = native._native_library()
        self.assertNotEqual(original, changed_headers)
        with mock.patch.object(native, '_NATIVE_FLAGS',
                               (*native._NATIVE_FLAGS, '-DNATIVE_TEST=1')):
            changed_flags = native._native_library()
        self.assertNotEqual(changed_headers, changed_flags)
        with mock.patch.object(native, '_compiler_command',
                               return_value=('/fake/clang++',)):
            changed_compiler = native._native_library()
        self.assertNotEqual(changed_headers, changed_compiler)
        self.assertEqual(self.run.call_count, 4)

    def test_tampered_binary_or_manifest_rebuilds(self):
        original = Path(native._native_library())
        original.write_bytes(b'not-ELF-any-more')
        self.assertEqual(native._native_library(), str(original))
        self.assertEqual(self.run.call_count, 2)
        original.with_suffix('.json').write_text('{invalid json')
        self.assertEqual(native._native_library(), str(original))
        self.assertEqual(self.run.call_count, 3)
        self.assertEqual(original.read_bytes(), b'\x7fELFmocked-native-library')

    def test_failed_compile_cleans_temporary_output_and_recovers(self):
        self.run.side_effect = lambda args, **kwargs: _fake_elf_run(
            args, fail=True, **kwargs)
        with self.assertRaises(subprocess.CalledProcessError):
            native._native_library()
        cache = self.root / 'target'
        self.assertEqual(list(cache.glob('*.so')), [])
        self.assertEqual(list(cache.glob('*.json')), [])
        self.assertEqual(list(cache.glob('.native-comparisons-*.so')), [])
        self.run.side_effect = _fake_elf_run
        self.assertTrue(Path(native._native_library()).is_file())

    def test_input_change_during_compile_refuses_stale_publication(self):
        def mutate(args, **kwargs):
            _fake_elf_run(args, **kwargs)
            self.source.write_text('int version() { return 99; }\n')
            return subprocess.CompletedProcess(args, 0)
        self.run.side_effect = mutate
        with self.assertRaisesRegex(RuntimeError, 'changed during compilation'):
            native._native_library()
        cache = self.root / 'target'
        self.assertFalse(list(cache.glob('*.so')))
        self.assertFalse(list(cache.glob('*.json')))

    @unittest.skipUnless(os.name == 'posix', 'requires Unix flock')
    def test_two_processes_compile_once_under_file_lock(self):
        # Local mock patches above do not propagate to a spawned process;
        # the worker installs its own isolated patches.
        context = multiprocessing.get_context('fork')
        event_file = self.root / 'compile-events.txt'
        queue = context.Queue()
        processes = [
            context.Process(target=_concurrent_caller,
                            args=(str(self.root), queue, str(event_file)))
            for _ in range(2)
        ]
        for process in processes:
            process.start()
        for process in processes:
            process.join(timeout=12)
            if process.is_alive():
                process.terminate()
                process.join(timeout=3)
            self.assertEqual(process.exitcode, 0)
        results = [queue.get(timeout=2) for _ in processes]
        self.assertTrue(all(kind == 'ok' for kind, _ in results), results)
        self.assertEqual(results[0][1], results[1][1])
        self.assertEqual(event_file.read_text().splitlines(), ['compile'])


if __name__ == '__main__':
    unittest.main()
