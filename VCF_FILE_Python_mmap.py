# 使用IO优化的策略彻底破产，考虑mmap(x)因为内存大小不同或者生成器
import gzip
import mmap

# 定义基类VCFFILE
class VCFFILE:
    def __init__(self,filename:str):

        # 私有变量
        self.__filename = filename
        self.__header_array = _set_header(filename)
        self.__genotype_array = _set_geno(filename)
    
def _set_header(filename):
    headers = []  # 用于存储所有以 "##" 开头的行
    with open(filename, "rt") as infile:  # 打开文件，自动处理缓冲
        for line in infile:
            if line.startswith("##"):
                headers.append(line)  # 将注释行添加到列表中
    return headers  # 返回收集的注释行


        



        
