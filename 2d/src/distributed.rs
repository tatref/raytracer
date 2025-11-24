use std::{
    fmt::Debug,
    io::{Read, Write},
    net::TcpStream,
};

use serde::{Deserialize, Serialize};

use crate::{Color, librt2d::World};

#[derive(Serialize, Deserialize)]
pub enum Message {
    Ping,
    Pong,
    JobRequest { world: World },
    JobResult { img: Vec<u8> },
    WorkRequest { columns: Vec<i32> },
    WorkResult { column: Vec<(i32, Color)> },
}

impl Debug for Message {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Message::Ping => write!(f, "Ping"),
            Message::Pong => write!(f, "Pong"),
            _ => unimplemented!(),
        }
    }
}

impl Message {
    pub fn send(&self, stream: &mut TcpStream) -> Result<usize, Box<dyn std::error::Error>> {
        let buf = rmp_serde::to_vec(&self)?;
        let size = (buf.len() as u64).to_le_bytes();

        stream.write_all(&size)?;
        stream.write_all(&buf)?;

        Ok(buf.len() + size.len())
    }

    pub fn recv(stream: &mut TcpStream) -> Result<Self, Box<dyn std::error::Error>> {
        let mut size = [0u8; 8];
        stream.read_exact(&mut size)?;
        let size = u64::from_le_bytes(size);

        let mut buf = vec![0u8; size as usize];
        stream.read_exact(&mut buf)?;

        let message: Message = rmp_serde::from_slice(&buf)?;
        Ok(message)
    }
}
