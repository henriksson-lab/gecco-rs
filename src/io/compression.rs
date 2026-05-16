//! Transparent decompression based on magic bytes.
//!
//! Counterpart of Python `gecco._meta.zopen`. The Python version supports
//! gzip, bzip2, xz, and lz4; the Rust port currently handles gzip and
//! returns the file as-is otherwise, which is enough for the formats GECCO
//! actually ships.

use std::io::{self, BufRead, BufReader, Read};

use flate2::read::GzDecoder;

const GZIP_MAGIC: &[u8] = &[0x1f, 0x8b];

/// Open a reader and transparently decompress it if its leading bytes
/// match a supported compression format.
///
/// Returns a boxed reader that yields the decompressed contents (or the
/// original bytes if the stream is not compressed).
pub fn zopen(reader: impl Read + 'static) -> io::Result<Box<dyn Read>> {
    let mut buf = BufReader::new(reader);
    let peek = buf.fill_buf()?;
    if peek.len() >= 2 && peek.starts_with(GZIP_MAGIC) {
        Ok(Box::new(GzDecoder::new(buf)))
    } else {
        Ok(Box::new(buf))
    }
}

/// Open a file at `path` and transparently decompress its contents,
/// mirroring the path-argument form of Python `gecco._meta.zopen`.
pub fn zopen_path(path: &std::path::Path) -> io::Result<Box<dyn Read>> {
    let file = std::fs::File::open(path)?;
    zopen(file)
}
