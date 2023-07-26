    # Conda软件安装目录，`conda env list`查看，如/anaconda3
    soft=/anaconda3/
    # 数据库database(db)位置，如管理员/db，个人~/db
    db=/db
    # 设置工作目录work directory(wd)，如meta
    wd=~/meta
    # 创建并进入工作目录
    mkdir -p $wd && cd $wd
    # 添加分析所需的软件、脚本至环境变量，添加至~/.bashrc中自动加载
    PATH=$PATH:$db/EasyMicrobiome/linux:$db/EasyMicrobiome/script