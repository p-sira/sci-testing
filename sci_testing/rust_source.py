"""Keep generated Rust reproducible under the installed formatter."""
import subprocess


def formatted(source):
    return subprocess.run(['rustfmt','--edition','2024','--emit','stdout'],input=source,
                          text=True,capture_output=True,check=True).stdout
