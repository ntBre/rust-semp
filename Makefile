BASE = /home/brent/Projects/rust-semp

eland:
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
	scp -C ${BASE}/target/release/rust-semp 'eland:semp/rust/.'

test:
	RUST_BACKTRACE=1 cargo test -- --test-threads=1 --nocapture

profile:
	RUSTFLAGS='-g' cargo build --release --bin profiling
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		${BASE}/target/release/profiling
