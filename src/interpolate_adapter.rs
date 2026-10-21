// Generated from the pinned Sci interpolation fixture dispatcher.
#![allow(dead_code, unused_imports)]
macro_rules! close {
    ($actual:expr,$expected:expr,$rtol:expr,$atol:expr $(,)?) => {{ emit($actual) }};
}
use crate::numeric::Numbers;
use std::cell::RefCell;
thread_local! { static OUTPUT: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) }; }
fn emit<T: Numbers>(value: T) {
    OUTPUT.with(|o| o.borrow_mut().extend(value.numbers().unwrap()));
}
use crate::common;
use common::close;
use common::{matrix, vector};

use sci::{Error, Spline, SurfaceSpline, interpolate as p};

fn spline(s: &str) -> Spline {
    let f: Vec<_> = s.split('|').collect();
    Spline {
        knots: vector(f[0]),
        coefficients: vector(f[1]),
        degree: f[2].parse().unwrap(),
    }
}
fn surface(s: &str) -> SurfaceSpline {
    let f: Vec<_> = s.split('|').collect();
    SurfaceSpline {
        x_knots: vector(f[0]),
        y_knots: vector(f[1]),
        coefficients: vector(f[2]),
        x_degree: f[3].parse().unwrap(),
        y_degree: f[4].parse().unwrap(),
    }
}
fn optional(s: &str) -> Option<Vec<f64>> {
    (!s.is_empty()).then(|| vector(s))
}
fn check(actual: &[f64], _expected: &str) {
    emit(actual);
}
fn check_matrix(actual: &[Vec<f64>], _expected: &str) {
    emit(actual);
}
fn extension(i: &str) -> p::SplineExtension {
    match i {
        "0" => p::SplineExtension::Extrapolate,
        "1" => p::SplineExtension::Zero,
        "2" => p::SplineExtension::Error,
        "3" => p::SplineExtension::Clamp,
        _ => panic!(),
    }
}
fn check_fit(s: &Spline, query: &[f64], expected: &[&str]) {
    check(
        &p::splev(query, s, 0, p::SplineExtension::Extrapolate).unwrap(),
        expected[1],
    );
    check(&s.knots, expected[2]);
    emit(s.degree);
}
fn check_transform(s: &Spline, query: &[f64], expected: &[&str]) {
    check(
        &p::splev(query, s, 0, p::SplineExtension::Extrapolate).unwrap(),
        expected[1],
    );
    check(&s.knots, expected[2]);
    check(&s.coefficients, expected[3]);
}

pub fn run(c: &[&str]) -> Vec<f64> {
    OUTPUT.with(|o| o.borrow_mut().clear());
    let a = &c[2..];
    let e = [""; 4];
    match c[1] {
        "make_interp_spline" => check_fit(
            &p::make_interp_spline(&vector(a[0]), &vector(a[1]), a[2].parse().unwrap()).unwrap(),
            &vector(a[3]),
            &e,
        ),
        "make_lsq_spline" => check_fit(
            &p::make_lsq_spline(
                &vector(a[0]),
                &vector(a[1]),
                &vector(a[2]),
                a[3].parse().unwrap(),
                optional(a[4]).as_deref(),
            )
            .unwrap(),
            &vector(a[5]),
            &e,
        ),
        "make_splrep" | "splrep" => {
            let f = if c[1] == "splrep" {
                p::splrep
            } else {
                p::make_splrep
            };
            check_fit(
                &f(
                    &vector(a[0]),
                    &vector(a[1]),
                    a[2].parse().unwrap(),
                    a[3].parse().unwrap(),
                    optional(a[4]).as_deref(),
                )
                .unwrap(),
                &vector(a[5]),
                &e,
            );
        }
        "make_smoothing_spline" => check_fit(
            &p::make_smoothing_spline(
                &vector(a[0]),
                &vector(a[1]),
                optional(a[2]).as_deref(),
                a[3].parse().unwrap(),
            )
            .unwrap(),
            &vector(a[4]),
            &e,
        ),
        "make_splprep" | "splprep" => {
            let f = if c[1] == "splprep" {
                p::splprep
            } else {
                p::make_splprep
            };
            let s = f(
                &matrix(a[0]),
                optional(a[1]).as_deref(),
                a[2].parse().unwrap(),
                a[3].parse().unwrap(),
                optional(a[4]).as_deref(),
            )
            .unwrap();
            let query = vector(a[5]);
            let mut out = vec![Vec::new(); query.len()];
            for component in &s.components {
                let v = p::splev(&query, component, 0, p::SplineExtension::Extrapolate).unwrap();
                for (row, value) in out.iter_mut().zip(v) {
                    row.push(value);
                }
            }
            check_matrix(&out, e[1]);
            check(&s.parameters, e[2]);
            check(&s.components[0].knots, e[3]);
            assert!(
                s.components
                    .iter()
                    .all(|c| c.knots == s.components[0].knots)
            );
        }
        "approximate_taylor_polynomial" => {
            let f: fn(f64) -> f64 = match a[0] {
                "square" => |x| x * x,
                "exp" => f64::exp,
                "sin" => f64::sin,
                "constant" => |_| 1.,
                _ => panic!(),
            };
            check(
                &p::approximate_taylor_polynomial(
                    f,
                    a[1].parse().unwrap(),
                    a[2].parse().unwrap(),
                    a[3].parse().unwrap(),
                    (!a[4].is_empty()).then(|| a[4].parse().unwrap()),
                )
                .unwrap(),
                e[1],
            );
        }
        "generate_knots" => check_matrix(
            &p::generate_knots(
                &vector(a[0]),
                &vector(a[1]),
                a[2].parse().unwrap(),
                a[3].parse().unwrap(),
                optional(a[4]).as_deref(),
            )
            .unwrap(),
            e[1],
        ),
        "splev" => check(
            &p::splev(
                &vector(a[0]),
                &spline(a[1]),
                a[2].parse().unwrap(),
                extension(a[3]),
            )
            .unwrap(),
            e[1],
        ),
        "splder" => check_transform(
            &p::splder(&spline(a[0]), a[1].parse().unwrap()).unwrap(),
            &vector(a[2]),
            &e,
        ),
        "splantider" => check_transform(
            &p::splantider(&spline(a[0]), a[1].parse().unwrap()).unwrap(),
            &vector(a[2]),
            &e,
        ),
        "spalde" => check_matrix(&p::spalde(&vector(a[0]), &spline(a[1])).unwrap(), e[1]),
        "splint" => close!(
            p::splint(a[0].parse().unwrap(), a[1].parse().unwrap(), &spline(a[2])).unwrap(),
            e[1].parse().unwrap(),
            2e-10,
            2e-10,
        ),
        "insert" => check_transform(
            &p::insert(a[0].parse().unwrap(), &spline(a[1]), a[2].parse().unwrap()).unwrap(),
            &vector(a[3]),
            &e,
        ),
        "sproot" => check(&p::sproot(&spline(a[0])).unwrap(), e[1]),
        "interpn" => check(
            &p::interpn(
                &matrix(a[0]),
                &vector(a[1]),
                &matrix(a[2]),
                if a[3] == "linear" {
                    p::GridMethod::Linear
                } else {
                    p::GridMethod::Nearest
                },
                a[4].parse().unwrap(),
                (!a[5].is_empty()).then(|| a[5].parse().unwrap()),
            )
            .unwrap(),
            e[1],
        ),
        "griddata" => check(
            &p::griddata(
                &matrix(a[0]),
                &vector(a[1]),
                &matrix(a[2]),
                match a[3] {
                    "linear" => p::ScatteredMethod::Linear,
                    "nearest" => p::ScatteredMethod::Nearest,
                    "cubic" => p::ScatteredMethod::Cubic,
                    _ => panic!(),
                },
                a[4].parse().unwrap(),
                a[5].parse().unwrap(),
            )
            .unwrap(),
            e[1],
        ),
        "bisplrep" => {
            let s = p::bisplrep(
                &vector(a[0]),
                &vector(a[1]),
                &vector(a[2]),
                &vector(a[3]),
                &vector(a[4]),
                a[5].parse().unwrap(),
                a[6].parse().unwrap(),
                optional(a[7]).as_deref(),
            )
            .unwrap();
            check(&s.x_knots, e[1]);
            check(&s.y_knots, e[2]);
            check(&s.coefficients, e[3]);
        }
        "bisplev" => check_matrix(
            &p::bisplev(
                &vector(a[0]),
                &vector(a[1]),
                &surface(a[2]),
                a[3].parse().unwrap(),
                a[4].parse().unwrap(),
            )
            .unwrap(),
            e[1],
        ),
        _ => panic!("unknown API {}", c[1]),
    };
    OUTPUT.with(|o| o.borrow().clone())
}
