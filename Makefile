eland:
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
	scp -C target/release/rust-semp 'eland:'

test:
	cargo test -- --test-threads=1 --nocapture
