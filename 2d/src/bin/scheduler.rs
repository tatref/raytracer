use std::{
    io::Write,
    net::{SocketAddr, TcpListener, TcpStream},
    sync::{
        Arc, Mutex,
        mpsc::{Receiver, Sender, SyncSender},
    },
    thread::{self, JoinHandle},
    time::Duration,
};

use raytracer::distributed::Message;

fn handle_client(mut stream: TcpStream, rx: Receiver<Message>, tx: SyncSender<Message>) {
    loop {
        thread::sleep(Duration::from_millis(100));

        // received message from network
        match Message::recv(&mut stream) {
            Ok(received) => println!("received message: {:?}", received),
            Err(e) => {
                let ioe = e.downcast::<std::io::Error>().unwrap();

                use std::io::ErrorKind::*;
                match ioe.kind() {
                    WouldBlock => {
                        // no data yet
                        println!("no data");
                    }
                    ConnectionAborted | ConnectionReset | UnexpectedEof => {
                        // client disconnected
                        println!("client disconnected: {:?}", ioe);
                    }
                    ioe => println!("unknown error: {:?}", ioe),
                }
            }
        }

        // receive message from main thread
        let Ok(message) = rx.try_recv() else {
            continue;
        };

        // send message to network
        message.send(&mut stream).unwrap();

        // send message to main thread
        // TODO
    }
}

fn main() {
    let socket: SocketAddr = "0.0.0.0:5000".parse().unwrap();
    let listener = TcpListener::bind(&socket).expect(&format!("Can't bind to socket {}", socket));

    let mut workers = Vec::new();

    loop {
        let (mut stream, worker_addr) = listener.accept().unwrap();
        stream.set_nonblocking(true).unwrap();

        let (tx1, rx1) = std::sync::mpsc::sync_channel::<Message>(1);
        let (tx2, rx2) = std::sync::mpsc::sync_channel::<Message>(1);
        println!("got worker {:?}", worker_addr);

        let t = thread::spawn(move || {
            handle_client(stream, rx1, tx2);
        });
        workers.push((worker_addr, t, tx1, rx2));
    }
}
