

module load compiler/cmake/3.24.1
module load compiler/devtoolset/7.3.1

module load compiler/rocm/dtk/25.04.3
module load mpi/hpcx/2.11.0/gcc-7.3.1

export HIP_PATH=/public/software/compiler/dtk/dtk-25.04.3
export ROCM_PATH=/public/software/compiler/dtk/dtk-25.04.3
export DEVICE_LIB_PATH=/public/software/compiler/dtk/dtk-25.04.3/lib/bitcode

# export CMAKE_PREFIX_PATH=/public/software/compiler/dtk/dtk-25.04.3:/public/software/compiler/dtk/dtk-25.04.3/lib64/cmake/amd_comgr:$CMAKE_PREFIX_PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/compiler/dtk/dtk-25.04.3/lib64
# export LIBRARY_PATH=/work1/baihe/users/sumeng/kmc/vendor/pkg/github.com/jbeder/yaml-cpp/lib:$LIBRARY_PATH
# export CPLUS_INCLUDE_PATH=/work1/baihe/users/sumeng/kmc/vendor/pkg/github.com/jbeder/yaml-cpp/include:$CPLUS_INCLUDE_PATH

export LIBRARY_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp/lib:$LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp/include:$CPLUS_INCLUDE_PATH
export CMAKE_PREFIX_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp:$CMAKE_PREFIX_PATH

# export LIBRARY_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp/lib:$LIBRARY_PATH
# export LIBRARY_PATH=$PWD/vendor/pkg/github.com/fmtlib/fmt/lib64:$LIBRARY_PATH
# export LIBRARY_PATH=$PWD/vendor/pkg/git.hpcer.dev/genshen/kiwi/lib/:$LIBRARY_PATH

# export CPLUS_INCLUDE_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp/include:$CPLUS_INCLUDE_PATH
# export CMAKE_PREFIX_PATH=$PWD/vendor/pkg/github.com/jbeder/yaml-cpp:$CMAKE_PREFIX_PATH

# export CPLUS_INCLUDE_PATH=$PWD/vendor/pkg/github.com/fmtlib/fmt/include:$CPLUS_INCLUDE_PATH
# export CPLUS_INCLUDE_PATH=$PWD/vendor/pkg/git.hpcer.dev/genshen/kiwi/include:$CPLUS_INCLUDE_PATH

# export CXXFLAGS="--gcc-toolchain=/opt/rh/devtoolset-7/root/usr"
# export LDFLAGS="--gcc-toolchain=/opt/rh/devtoolset-7/root/usr"