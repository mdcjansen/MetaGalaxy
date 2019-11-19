FROM ubuntu:18.04 AS run-env

LABEL Maintainer="Rob Egan<RSEgan@lbl.gov>"

WORKDIR /root

# This is necessary because the upgrade sometimes prompts for input
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update  && \
    apt-get install -y libgomp1 && \
    apt-get autoremove -y && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/*

FROM run-env as builder

RUN apt-get update && \
    apt-get install -y build-essential autoconf libboost-all-dev cmake git curl libncurses5-dev zlib1g-dev

# copy the git tree (minus Dockerfile) to metabat subdir
COPY . metabat

RUN cd metabat && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && \
    make -j8 && \
    make install && \
    cd .. && \
    rm -rf build

FROM run-env
WORKDIR /root
ENV PATH=$PATH:/root/bin
COPY --from=builder /usr/local /usr/local

env PATH=/usr/local/bin:$PATH

CMD ["/usr/local/bin/runMetaBat.sh"]

# build and deploy with this command
# docker build --tag robegan21/metabat:$(git describe --tags) --tag robegan21/metabat:latest . && docker push robegan21/metabat
