use std::net::{SocketAddr, TcpStream};

fn main() {
    let args: Vec<_> = std::env::args().collect();

    if args.len() != 2 {
        println!("Usage: {} <scheduler_addr>", args[0]);
        panic!();
    }

    let scheduler_addr: SocketAddr = args[1].parse().unwrap();
    let mut socket =
        TcpStream::connect_timeout(&scheduler_addr, std::time::Duration::from_secs(2)).unwrap();

    loop {}
}
