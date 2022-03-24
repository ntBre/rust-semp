BASE = /home/brent/Projects/rust-semp
ARGS =
TESTFLAGS = --test-threads=1 --nocapture
SHORT = 0

ifeq ($(SHORT),0)
TESTFLAGS += --include-ignored
endif

eland:
	RUSTFLAGS='-C target-feature=+crt-static' cargo build --release
	scp -C ${BASE}/target/release/rust-semp 'eland:semp/rust/.'

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
