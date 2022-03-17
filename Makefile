eland:
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
	scp -C target/release/rust-semp 'eland:'
