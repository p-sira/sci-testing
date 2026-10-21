//! Shared typed argument conversion and numeric output ordering for registrations.
use super::{Result, scalar};
use sci::complex::Complex;
use serde_json::Value;

pub trait Argument: Sized {
    fn read(value: &Value) -> Result<Self>;
}
impl Argument for f64 {
    fn read(value: &Value) -> Result<Self> {
        scalar(value)
    }
}
macro_rules! integers {
    ($($kind:ty),*) => {$(impl Argument for $kind {
        fn read(value: &Value) -> Result<Self> {
            if let Some(number) = value.as_u64() {
                return <$kind>::try_from(number).map_err(|_| "integer out of range".into());
            }
            if let Some(number) = value.as_i64() {
                return <$kind>::try_from(number).map_err(|_| "integer out of range".into());
            }
            let number = scalar(value)?;
            if !number.is_finite() || number.fract() != 0.0 {
                return Err("expected a finite integer".into());
            }
            <$kind>::try_from(number as i128).map_err(|_| "integer out of range".into())
        }
    })*};
}
integers!(u32, i32, usize, u64, i64, isize);
impl Argument for bool {
    fn read(value: &Value) -> Result<Self> {
        value.as_bool().ok_or_else(|| "expected bool".into())
    }
}
impl Argument for Complex {
    fn read(value: &Value) -> Result<Self> {
        Ok(Self {
            re: scalar(&value[0])?,
            im: scalar(&value[1])?,
        })
    }
}

pub trait Numbers {
    fn numbers(self) -> Result<Vec<f64>>;
}
macro_rules! scalars {
    ($($kind:ty),*) => {$(impl Numbers for $kind {
        fn numbers(self) -> Result<Vec<f64>> { Ok(vec![self as f64]) }
    })*};
}
scalars!(f64, usize, u32, i32, u64, i64);
impl Numbers for Complex {
    fn numbers(self) -> Result<Vec<f64>> {
        Ok(vec![self.re, self.im])
    }
}
impl<T: Numbers> Numbers for Vec<T> {
    fn numbers(self) -> Result<Vec<f64>> {
        self.into_iter()
            .map(Numbers::numbers)
            .collect::<Result<Vec<_>>>()
            .map(|rows| rows.concat())
    }
}
impl<T: Numbers, const N: usize> Numbers for [T; N] {
    fn numbers(self) -> Result<Vec<f64>> {
        Vec::from(self).numbers()
    }
}
impl<T: Numbers, E: std::fmt::Display> Numbers for std::result::Result<T, E> {
    fn numbers(self) -> Result<Vec<f64>> {
        self.map_err(|error| error.to_string())?.numbers()
    }
}
macro_rules! tuples {
    ($($kind:ident:$index:tt),*) => {
        impl<$($kind: Numbers),*> Numbers for ($($kind,)*) {
            fn numbers(self) -> Result<Vec<f64>> { Ok([$((self.$index).numbers()?,)*].concat()) }
        }
    };
}
tuples!(A:0, B:1);
tuples!(A:0, B:1, C:2, D:3);
tuples!(A:0, B:1, C:2, D:3, E:4, F:5, G:6, H:7);

impl Argument for sci::special::ErrorPolicy {
    fn read(value: &Value) -> Result<Self> {
        match scalar(value)? {
            0.0 => Ok(Self::Ieee),
            1.0 => Ok(Self::Finite),
            _ => Err("unknown special error policy".into()),
        }
    }
}

impl<T: Argument> Argument for Vec<T> {
    fn read(value: &Value) -> Result<Self> {
        value
            .as_array()
            .ok_or("expected array")?
            .iter()
            .map(T::read)
            .collect()
    }
}
impl<T: Argument> Argument for Option<T> {
    fn read(value: &Value) -> Result<Self> {
        if value.is_null() {
            Ok(None)
        } else {
            T::read(value).map(Some)
        }
    }
}

impl<T: Numbers + Clone> Numbers for &[T] {
    fn numbers(self) -> Result<Vec<f64>> {
        self.to_vec().numbers()
    }
}
impl<T: Numbers + Clone> Numbers for &Vec<T> {
    fn numbers(self) -> Result<Vec<f64>> {
        self.as_slice().numbers()
    }
}
impl Numbers for bool {
    fn numbers(self) -> Result<Vec<f64>> {
        Ok(vec![if self { 1. } else { 0. }])
    }
}
impl Numbers for String {
    fn numbers(self) -> Result<Vec<f64>> {
        Ok(std::iter::once(self.len() as f64)
            .chain(self.bytes().map(f64::from))
            .collect())
    }
}
