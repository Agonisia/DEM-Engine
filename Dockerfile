FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# 更新包列表并安装所有依赖
RUN apt-get update && apt-get install -y \
    # 构建工具
    build-essential \
    cmake \
    git \
    wget \
    # Python和科学计算包
    python3 \
    python3-dev \
    python3-numpy \
    python3-matplotlib \
    python3-scipy \
    python3-pandas \
    python3-h5py \
    # MPI支持
    libopenmpi-dev \
    openmpi-bin \
    # 数学库
    libeigen3-dev \
    libboost-all-dev \
    liblapack-dev \
    libblas-dev \
    libgsl-dev \
    # 其他工具
    vim \
    htop \
    && rm -rf /var/lib/apt/lists/*

# 设置工作目录
WORKDIR /workspace

# 设置环境变量
ENV OMP_NUM_THREADS=8
ENV PYTHONPATH=/workspace:$PYTHONPATH

# 创建结果目录
RUN mkdir -p /results

CMD ["/bin/bash"]