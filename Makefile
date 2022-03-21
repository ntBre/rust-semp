BASE = /home/brent/Projects/rust-semp

eland:
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
	scp -C ${BASE}/target/release/rust-semp 'eland:semp/rust/.'

test:
	RUST_BACKTRACE=1 cargo test -- --test-threads=1 --nocapture
