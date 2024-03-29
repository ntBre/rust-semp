BASE = /home/brent/Projects/rust-semp
ARGS =
TESTFLAGS = -- --test-threads=1 --nocapture
SHORT = 0

ifeq ($(SHORT),0)
TESTFLAGS += --include-ignored
endif

ELAND_DEST = 'eland:programs/semp/.'
WOODS_DEST = 'woods:bin/semp'

clippy:
	cargo clippy --workspace --tests

build:
# see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
	RUSTFLAGS='-C target-feature=+crt-static' \
	cargo build --release --target x86_64-unknown-linux-gnu --bin rust-semp

eland: build
	scp -C ${BASE}/target/x86_64-unknown-linux-gnu/release/rust-semp ${ELAND_DEST}

woods: build
	scp -C ${BASE}/target/x86_64-unknown-linux-gnu/release/rust-semp ${WOODS_DEST}${ALPHA}

eland.scripts: scripts/time.awk
	scp -C $? ${ELAND_DEST}

test:
	RUST_BACKTRACE=1 cargo test ${TESTFLAGS} ${ARGS}

clean:
	cargo clean

profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
		--collect-jumps=yes --simulate-cache=yes		\
		${BASE}/target/release/$(1)
cover:
	cargo tarpaulin --skip-clean --timeout 3600 --out html \
	    --exclude-files src/bin/* src/main.rs \
            --color=never \
	    ${TESTFLAGS} ${ARGS}

profile.one_iter:
	$(call profile,one_iter)

profile.num_jac:
	$(call profile,num_jac)

profile.full:
	$(call profile,full)

profile.freq_num_jac:
	$(call profile,freq_num_jac)

memprofile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
		heaptrack ${BASE}/target/release/$(1)

memprofile.full:
	$(call memprofile,full)

memprofile.freq_num_jac:
	$(call memprofile,freq_num_jac)

memprofile.norm_num_jac:
	$(call memprofile,norm_num_jac)

time:
	RUSTFLAGS='-g' cargo build --release --bin full
	sh -c "time ${BASE}/target/release/full"
