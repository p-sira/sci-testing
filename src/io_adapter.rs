//! Public Sci I/O comparisons against serialized, independent SciPy fixtures.
use sci::io::{NumericData, matlab, matrix_market, wavfile};
use serde_json::Value;
use std::error::Error;
use std::io::Cursor;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn bytes(a: &Value) -> Result<Vec<u8>> {
    let hex = a["payload_hex"].as_str().ok_or("missing I/O hex input")?;
    if !hex.len().is_multiple_of(2) {
        return Err("I/O hex input has odd length".into());
    }
    (0..hex.len())
        .step_by(2)
        .map(|i| Ok(u8::from_str_radix(&hex[i..i + 2], 16)?))
        .collect()
}

// Every sample is (class, value, zero-sign). Classes are 0 for finite, +1/-1
// for infinity, and 2 for NaN. Integer64 payloads additionally split into
// lossless upper/lower 32-bit words, each carrying its own class tag.
fn scalar(out: &mut Vec<f64>, x: f64, signed_zero: bool) {
    if x.is_nan() {
        out.extend([2., 0., 0.]);
    } else if x.is_infinite() {
        out.extend([x.signum(), 0., 0.]);
    } else {
        out.extend([
            0.,
            x,
            if signed_zero && x == 0. && x.is_sign_negative() {
                -1.
            } else {
                1.
            },
        ]);
    }
}

fn numeric(out: &mut Vec<f64>, samples: &NumericData) {
    let code = match samples {
        NumericData::I8(_) => 1.,
        NumericData::U8(_) => 2.,
        NumericData::I16(_) => 3.,
        NumericData::U16(_) => 4.,
        NumericData::I32(_) => 5.,
        NumericData::U32(_) => 6.,
        NumericData::I64(_) => 7.,
        NumericData::U64(_) => 8.,
        NumericData::F32(_) => 9.,
        NumericData::F64(_) => 10.,
    };
    out.extend([code, samples.len() as f64]);
    match samples {
        NumericData::I64(v) => {
            for &value in v {
                let bits = value as u64;
                scalar(out, (bits >> 32) as f64, true);
                scalar(out, (bits & 0xffff_ffff) as f64, true);
            }
        }
        NumericData::U64(v) => {
            for &bits in v {
                scalar(out, (bits >> 32) as f64, true);
                scalar(out, (bits & 0xffff_ffff) as f64, true);
            }
        }
        _ => {
            for value in samples.to_f64() {
                scalar(out, value, true);
            }
        }
    }
}

fn text(out: &mut Vec<f64>, s: &str) {
    out.push(s.len() as f64);
    out.extend(s.bytes().map(f64::from));
}

fn market_info(info: matrix_market::MatrixMarketInfo) -> Vec<f64> {
    let format = match info.format {
        "array" => 0.,
        "coordinate" => 1.,
        _ => unreachable!(),
    };
    let field = match info.field {
        "real" => 0.,
        "integer" => 1.,
        "pattern" => 2.,
        _ => unreachable!(),
    };
    let symmetry = match info.symmetry {
        "general" => 0.,
        "symmetric" => 1.,
        "skew-symmetric" => 2.,
        "hermitian" => 3.,
        _ => unreachable!(),
    };
    vec![
        info.rows as f64,
        info.cols as f64,
        info.entries as f64,
        format,
        field,
        symmetry,
    ]
}

fn market_matrix(a: &matrix_market::MarketMatrix) -> Result<Vec<f64>> {
    let (is_sparse, rows, cols, values) = match a {
        matrix_market::MarketMatrix::Dense { rows, cols, data } => (0., *rows, *cols, data.clone()),
        matrix_market::MarketMatrix::Sparse(matrix) => {
            let values = sci::sparse::construct::to_dense(matrix)?;
            (
                1.,
                matrix.rows,
                matrix.cols,
                values.into_iter().flatten().collect(),
            )
        }
    };
    let mut out = vec![is_sparse, rows as f64, cols as f64];
    for value in values {
        // SciPy's Matrix Market text parser does not preserve textual -0.
        scalar(&mut out, value, false);
    }
    Ok(out)
}

fn wav_audio(audio: &wavfile::Audio) -> Vec<f64> {
    let mut out = vec![
        audio.sample_rate as f64,
        audio.channels as f64,
        audio.bits_per_sample as f64,
    ];
    numeric(&mut out, &audio.samples);
    out
}

fn mat_class(name: &str) -> Result<f64> {
    Ok(match name {
        "double" => 1.,
        "single" => 2.,
        "int8" => 3.,
        "uint8" => 4.,
        "int16" => 5.,
        "uint16" => 6.,
        "int32" => 7.,
        "uint32" => 8.,
        "int64" => 9.,
        "uint64" => 10.,
        _ => return Err(format!("unsupported MAT class {name}").into()),
    })
}

fn mat_metadata(info: &[matlab::MatVariableInfo]) -> Result<Vec<f64>> {
    let mut out = vec![info.len() as f64];
    for variable in info {
        text(&mut out, &variable.name);
        out.push(variable.shape.len() as f64);
        out.extend(variable.shape.iter().map(|&v| v as f64));
        out.extend([
            mat_class(variable.class)?,
            f64::from(u8::from(variable.complex)),
        ]);
    }
    Ok(out)
}

fn mat_values(file: &matlab::MatFile) -> Vec<f64> {
    let mut out = vec![file.len() as f64];
    for (name, array) in file {
        text(&mut out, name);
        out.push(array.shape.len() as f64);
        out.extend(array.shape.iter().map(|&v| v as f64));
        out.push(f64::from(u8::from(array.imag.is_some())));
        numeric(&mut out, &array.real);
        if let Some(imag) = &array.imag {
            numeric(&mut out, imag);
        }
    }
    out
}

fn write_bytes(api: &str, source: &[u8], args: &Value) -> Result<Vec<u8>> {
    match api {
        "io.matrix_market.mmwrite" => {
            let matrix = matrix_market::mmread(source)?;
            let mut out = Vec::new();
            matrix_market::mmwrite(&mut out, &matrix)?;
            Ok(out)
        }
        "io.wavfile.write" => {
            let audio = wavfile::read(source)?;
            let mut out = Cursor::new(Vec::new());
            wavfile::write(&mut out, &audio)?;
            Ok(out.into_inner())
        }
        "io.matlab.savemat" => {
            let variables = matlab::loadmat(source)?;
            let compressed = args["compressed"]
                .as_bool()
                .ok_or("missing MAT compression")?;
            let mut out = Vec::new();
            matlab::savemat(&mut out, &variables, compressed)?;
            Ok(out)
        }
        _ => Err(format!("I/O byte probe requires a writer: {api}").into()),
    }
}

fn evaluate_success(api: &str, data: &[u8], a: &Value) -> Result<Vec<f64>> {
    Ok(match api {
        "io.matrix_market.mminfo" => market_info(matrix_market::mminfo(data)?),
        "io.matrix_market.mmread" => market_matrix(&matrix_market::mmread(data)?)?,
        "io.matrix_market.mmwrite" => market_matrix(&matrix_market::mmread(
            write_bytes(api, data, a)?.as_slice(),
        )?)?,
        "io.wavfile.read" => wav_audio(&wavfile::read(data)?),
        "io.wavfile.write" => wav_audio(&wavfile::read(write_bytes(api, data, a)?.as_slice())?),
        "io.wavfile.signal" => {
            let audio = wavfile::read(data)?;
            let channel = a["channel"].as_u64().ok_or("missing WAV channel")? as usize;
            let signal = wavfile::signal(&audio, channel)?;
            let mut values = vec![signal.len() as f64];
            for sample in signal {
                scalar(&mut values, sample, true);
            }
            values
        }
        "io.matlab.matfile_version" => {
            let (major, minor) = matlab::matfile_version(data)?;
            vec![major as f64, minor as f64]
        }
        "io.matlab.whosmat" => mat_metadata(&matlab::whosmat(data)?)?,
        "io.matlab.loadmat" => mat_values(&matlab::loadmat(data)?),
        "io.matlab.savemat" => mat_values(&matlab::loadmat(write_bytes(api, data, a)?.as_slice())?),
        _ => return Err(format!("unregistered I/O API: {api}").into()),
    })
}

/// Runs a Sci public I/O API. Reject cases compare actual error behavior as 1/0.
pub fn evaluate(api: &str, args: &Value) -> Result<Vec<f64>> {
    let data = bytes(args)?;
    match args["operation"].as_str().ok_or("missing I/O operation")? {
        "success" => evaluate_success(api, &data, args),
        "reject" => Ok(vec![f64::from(u8::from(
            evaluate_success(api, &data, args).is_err(),
        ))]),
        // Test-only diagnostic: the normal accuracy path always uses success
        // and reject. Return actual Rust writer bytes as exact integer values
        // for a separate Python test to read with public SciPy codecs.
        "emit_bytes" => Ok(write_bytes(api, &data, args)?
            .into_iter()
            .map(f64::from)
            .collect()),
        operation => Err(format!("invalid I/O operation: {operation}").into()),
    }
}
