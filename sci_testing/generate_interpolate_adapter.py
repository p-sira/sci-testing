"""Reuse the pinned interpolation dispatch and record its observable outputs."""
from .rust_source import formatted
from .api_registry import SCI

from . import ROOT


def generated():
    text=(SCI/'tests/interpolation.rs').read_text()
    prefix=text[:text.index('#[test]')].replace('mod common;','use crate::common;')
    start=prefix.index('fn check(actual:');end=prefix.index('fn extension(',start)
    prefix=prefix[:start]+'''fn check(actual: &[f64], _expected: &str) { emit(actual); }
fn check_matrix(actual: &[Vec<f64>], _expected: &str) { emit(actual); }
'''+prefix[end:]
    prefix=prefix.replace('assert_eq!(s.degree, expected[3].parse::<usize>().unwrap());','emit(s.degree);')
    start=text.index('        let result = std::panic::catch_unwind(|| match c[1] {')
    end=text.index('        if result.is_err()',start)
    body=text[start:end].replace('close(', 'close!(').replace('let result = std::panic::catch_unwind(|| match c[1] {','match c[1] {').rsplit('});',1)[0]+'};'
    return '''// Generated from the pinned Sci interpolation fixture dispatcher.
#![allow(dead_code,unused_imports)]
macro_rules! close { ($actual:expr,$expected:expr,$rtol:expr,$atol:expr $(,)?) => {{ emit($actual) }}; }
use crate::numeric::Numbers;
use std::cell::RefCell;
thread_local! { static OUTPUT: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) }; }
fn emit<T: Numbers>(value: T) { OUTPUT.with(|o|o.borrow_mut().extend(value.numbers().unwrap())); }
'''+prefix+'''pub fn run(c: &[&str]) -> Vec<f64> {
    OUTPUT.with(|o|o.borrow_mut().clear());
    let a=&c[2..];
    let e=["";4];
'''+body+'''
    OUTPUT.with(|o|o.borrow().clone())
}
'''

if __name__=='__main__':(ROOT/'src/interpolate_adapter.rs').write_text(formatted(generated()))
