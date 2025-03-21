From ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive

# RUN apt-get --allow-releaseinfo-change update
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils 
RUN apt-get install -f --fix-broken -y --no-install-recommends \
    locales \
    python3-pip \
    ca-certificates \
    netbase \
    gnupg \
    dirmngr \
    git \
    openssh-client \
    rpm \
    sshpass \
    valgrind \
    gdb \
    clang-format \
    ninja-build \
    doxygen \
    graphviz \
    curl \
    wget \
    build-essential \
    cmake \
    libboost-all-dev \
    pkg-config \
    make \
    cmake \
    g++ \
    libeigen3-dev \
    gfortran \
    libmetis-dev \
    libopenvdb-dev \
    libboost-all-dev \
    libtbb2 \
    libtbb-dev \
    libz-dev \
    libgmp3-dev \
    libmpc-dev \
    vim \
    libeigen3-dev \
    libcgal-dev \
    libblas-dev liblapack-dev \
    xorg-dev libglu1-mesa-dev \
    libglu1-mesa-dev freeglut3-dev mesa-common-dev \
    libgmsh-dev gmsh\
    libfmt-dev libglm-dev libopenblas-dev liblapack-dev liblapack3 libopenblas-base \
    && apt-get clean

RUN apt-get update && apt-get upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates gpg-agent software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install mesa-utils -y --no-install-recommends  && rm -rf /var/lib/apt/lists/*
# repository to install Intel(R) oneAPI Libraries
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | gpg --dearmor | tee /usr/share/keyrings/intel-oneapi-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/intel-oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main " > /etc/apt/sources.list.d/oneAPI.list

# repository to install Intel(R) GPU drivers
RUN curl -fsSL https://repositories.intel.com/graphics/intel-graphics.key | gpg --dearmor | tee /usr/share/keyrings/intel-graphics-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/intel-graphics-archive-keyring.gpg arch=amd64] https://repositories.intel.com/graphics/ubuntu jammy flex" > /etc/apt/sources.list.d/intel-graphics.list

RUN apt-get update && apt-get upgrade -y && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  ca-certificates build-essential pkg-config gnupg libarchive13 openssh-server openssh-client wget net-tools git intel-oneapi-mkl-devel-2024.2 intel-oneapi-tbb-devel-2022.0 && \
  rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/liyuesolo/SuiteSparse && cd SuiteSparse && git checkout dev && cd build && cmake .. && make -j8 install && cd ../..

RUN cp -r /usr/include/eigen3/* /usr/local/include/

ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES all

# RUN git clone --recursive https://github.com/coin-or/CppAD.git && cd CppAD && mkdir build && cd build && cmake .. -G Ninja && cmake --build . -j8 && ninja install -j8 && cd ../..

# RUN apt-get update && apt-get install libgtest-dev -y && mkdir -p $HOME/build && cd $HOME/build &&  cmake /usr/src/googletest/googletest && make -j8 &&  cp lib/libgtest* /usr/lib/ && cd .. &&  rm -rf build &&  mkdir /usr/local/lib/googletest &&  ln -s /usr/lib/libgtest.a /usr/local/lib/googletest/libgtest.a &&  ln -s /usr/lib/libgtest_main.a /usr/local/lib/googletest/libgtest_main.a
# RUN git clone --recursive https://github.com/joaoleal/CppADCodeGen.git && cd CppADCodeGen && mkdir build && cd build && cmake .. -G Ninja && cmake --build . -j8 && ninja install  -j8 && cd ../..

# For NeuralMaterialModel
# RUN wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-2.15.0.tar.gz \
#   && tar -C /usr/local -xzf libtensorflow-cpu-linux-x86_64-2.15.0.tar.gz
# RUN git clone --recurse-submodules https://github.com/serizba/cppflow && cd cppflow && mkdir build && cd build && cmake .. && make -j8 install && cd ../..
