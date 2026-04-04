//! Transparent decompression based on magic bytes.

use std::io::{self, BufRead, BufReader, Read};

use flate2::read::GzDecoder;

const GZIP_MAGIC: &[u8] = &[0x1f, 0x8b];

/// Open a reader and transparently decompress if gzipped.
/// Returns a boxed reader.
pub fn zopen(reader: impl Read + 'static) -> io::Result<Box<dyn Read>> {
    let mut buf = BufReader::new(reader);
    let peek = buf.fill_buf()?;
    if peek.len() >= 2 && peek.starts_with(GZIP_MAGIC) {
        Ok(Box::new(GzDecoder::new(buf)))
    } else {
        Ok(Box::new(buf))
    }
}

/// Open a file and transparently decompress.
pub fn zopen_path(path: &std::path::Path) -> io::Result<Box<dyn Read>> {
    let file = std::fs::File::open(path)?;
    zopen(file)
}
