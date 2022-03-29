BASE = /home/brent/Projects/rust-semp
ARGS =
TESTFLAGS = --test-threads=1 --nocapture
SHORT = 0

ifeq ($(SHORT),0)
TESTFLAGS += --include-ignored
endif

ELAND_DEST = 'eland:programs/semp/.'
eland:
# see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release --target x86_64-unknown-linux-gnu
	scp -C ${BASE}/target/x86_64-unknown-linux-gnu/release/rust-semp ${ELAND_DEST}

test:
	RUST_BACKTRACE=1 cargo test -- ${TESTFLAGS} ${ARGS}

profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		${BASE}/target/release/$(1)

profile.one_iter:
	$(call profile,one_iter)

profile.num_jac:
	$(call profile,num_jac)

profile.full:
	$(call profile,full)
