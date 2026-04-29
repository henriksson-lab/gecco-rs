use std::fmt;
use std::io::{self, Write};
use std::sync::{Mutex, PoisonError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OutputStream {
    Stdout,
    Stderr,
}

pub trait RunOutput: Send + Sync {
    fn write(&self, stream: OutputStream, args: fmt::Arguments<'_>) -> io::Result<()>;

    fn stdout(&self, args: fmt::Arguments<'_>) -> io::Result<()> {
        self.write(OutputStream::Stdout, args)
    }

    fn stderr(&self, args: fmt::Arguments<'_>) -> io::Result<()> {
        self.write(OutputStream::Stderr, args)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct StdioOutput;

impl RunOutput for StdioOutput {
    fn write(&self, stream: OutputStream, args: fmt::Arguments<'_>) -> io::Result<()> {
        match stream {
            OutputStream::Stdout => {
                let mut out = io::stdout().lock();
                out.write_fmt(args)?;
                out.write_all(b"\n")
            }
            OutputStream::Stderr => {
                let mut err = io::stderr().lock();
                err.write_fmt(args)?;
                err.write_all(b"\n")
            }
        }
    }
}

pub struct SharedWriterOutput<W> {
    writer: Mutex<W>,
    label_streams: bool,
}

impl<W> SharedWriterOutput<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer: Mutex::new(writer),
            label_streams: false,
        }
    }

    pub fn with_stream_labels(writer: W) -> Self {
        Self {
            writer: Mutex::new(writer),
            label_streams: true,
        }
    }

    pub fn into_inner(self) -> Result<W, PoisonError<W>> {
        self.writer.into_inner()
    }
}

impl<W: Write + Send> RunOutput for SharedWriterOutput<W> {
    fn write(&self, stream: OutputStream, args: fmt::Arguments<'_>) -> io::Result<()> {
        let mut writer = self
            .writer
            .lock()
            .map_err(|_| io::Error::other("run output writer lock poisoned"))?;
        if self.label_streams {
            match stream {
                OutputStream::Stdout => writer.write_all(b"[stdout] ")?,
                OutputStream::Stderr => writer.write_all(b"[stderr] ")?,
            }
        }
        writer.write_fmt(args)?;
        writer.write_all(b"\n")
    }
}
