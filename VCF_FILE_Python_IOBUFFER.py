import gzip
from io import BufferedReader, BufferedWriter

class VCFFILE:
    def __init__(self, filename: str):
        self.__filename = filename
        self.__vcftext = self._set_vcftext()

    def _set_vcftext(self):
        header_count = 0
        vcf_lines = []
        def _is_gzipped() -> bool:
            with open(self.__filename, "rb") as in_binfile:
                magic = in_binfile.read(2)
                return magic == b'\x1f\x8b'
        if _is_gzipped():
            with gzip.open(self.__filename, "rt", buffering=4096) as in_file:
                reader = BufferedReader(in_file)
                for line in reader:
                    if line.startswith('##'):
                        header_count += 1
                    elif line.startswith("#CHROM"):
                        vcf_lines.append(line)
                        break
                    else:
                        vcf_lines.append(line)
        else:
            with open(self.__filename, "rt", buffering=4096) as in_file:
                reader = BufferedReader(in_file)
                for line in reader:
                    if line.startswith('##'):
                        header_count += 1
                    elif line.startswith("#CHROM"):
                        vcf_lines.append(line)
                        break
                    else:
                        vcf_lines.append(line)
        return vcf_lines

    def vcf(self):
        return self.__vcftext

# 对于写入操作，使用BufferedWriter
def write_to_file(data, filename):
    with gzip.open(filename, "wt", buffering=4096) as out_file:
        writer = BufferedWriter(out_file)
        for line in data:
            writer.write(line)
            writer.flush()  # 确保数据被写入

# 定义现代人VCF子类

class MODERNVCF(VCFFILE):
    def __init__(self,filename:str):
        super().__init__(filename)
        # 私有变量
        self.__index_dict:dict = self._set_index()

        # 公有变量
        self.pos_pool:set = self._set_pool()


    # 私有变量初始化函数：
    def _set_index(self):
        header = self.vcf()[0].strip().split("\t")
        index_out = {}
        try:
            index_out.update({
                "chrom":header.index("#CHROM"),
                "pos":header.index("POS"),
                "ref":header.index("REF"),
                "alt":header.index("ALT"),
                "sample":header.index("FORMAT")+1
            })
        except ValueError:
            print("Header Error, Cannot find right index,Check your header in Modern VCF!\n")
            exit(1)        
        return index_out        

    # 公有变量初始化函数
    def _set_pool(self):
        pos_pool = set()
        pos = self.__index_dict["pos"]
        for line in self.vcf()[1:]:
            line = line.strip().split("\t")
            pos_pool.add(line[pos])
        return pos_pool
        

    # 私有变量访问接口
    def get(self,attr:str):
        if attr in self.__index_dict.keys():
            return self.__index_dict[attr]
        else:
            print(f"The Porperty {attr} is not present in the modern VCF file,Please Check!\n")
            exit(1)
        
# 定义古人VCF子类
class ARCHAICVCF(VCFFILE):
    def __init__(self,filename:str):
        super().__init__(filename)
        # 私有变量
        self.__index_dict = self._set_index()

    
    # 私有变量初始化函数
    def _set_index(self):
        index_out = {}
        header = self.vcf()[0].strip().split("\t")
        try:
            index_out.update({
                "chrom":header.index("#CHROM"),
                "pos":header.index("POS"),
                "ref":header.index("REF"),
                "alt":header.index("ALT"),
                "archaic":header.index("FORMAT")+1
            })
        except ValueError:
            print("Header Error, Cannot find right index,Check your header in Archaic VCF!\n")
            exit(1)        
        return index_out

    # 私有变量访问接口
    def get(self,attr:str):
        if attr in self.__index_dict.keys():
            return self.__index_dict[attr]
        else:
            print(f"The Porperty {attr} is not present in the archaic VCF file,Please Check!\n")
            exit(1)     

    # 重写+号操作符用于merge操作，并返回状态
    def __add__(self,modern_file:MODERNVCF) -> bool:
        header = self.vcf()[0].strip().split("\t")
        data_sample = self.vcf()[1].strip().split("\t")
        data_modern = modern_file.vcf()[1].strip().split("\t")
        archaic = header[self.get("archaic")]
        chrom = data_sample[self.get("chrom")]
        chrom_modern = data_modern[self.get("chrom")]
        if not chrom_modern == chrom:
            print("You input the Worng files, chromsomes of Archaic VCF and Modern VCF cannot match!\n")
            exit(1)
        sample = modern_file.vcf()[0].strip().split("\t")[modern_file.get("sample"):]
        # 输出header
        filename = f"{archaic}_{chrom}.gz"
        write_to_file("\t".join(['chrom', 'pos', 'ref', 'alt', archaic] + sample )+"\n",filename)
        # 初始化现代人VCF的文件指针
        seek_pointer = 1
        # 遍历古人中的行
        for line in self.vcf()[1:]:
            line = line.strip().split("\t")
            pos = int(line[self.get("pos")])
            ref = line[self.get("ref")]
            alt = line[self.get("alt")]
            if not (len(ref) == 1 and len(alt) ==1):
                continue
            # 转换古人的基因型
            try:
                archaic_geno = str(int(line[self.get("archaic")][0])+int(line[self.get("archaic")][2]))
            except ValueError:
                archaic_geno = "0"
            # 检查是否在archaic_genotype中非0，modern中不出现的
            if not archaic_geno == "0" and pos not in modern_file.pos_pool:
                write_to_file("\t".join([chrom,str(pos),ref,alt,archaic_geno]+["0"]*len(sample))+"\n",filename)
                continue
            
            # 在现代人vcf文件中遍历寻找
            if seek_pointer <len(modern_file.vcf())-1:
                print("Right!\n")
                for md_line in modern_file.vcf()[seek_pointer:]:
                    md_line = md_line.strip().split("\t")
                    md_pos = int(md_line[modern_file.get("pos")])
                    md_ref = md_line[modern_file.get("ref")]
                    md_alt = md_line[modern_file.get("alt")]
                    if pos > md_pos:
                        seek_pointer += 1
                        continue
                    elif pos == md_pos:
                        if (len(md_ref) ==1 and len(md_alt) ==1 and md_ref == ref 
                            and (md_alt ==alt or alt == ".")):
                            md_geno = [str(int(x[0])+int(x[2])) 
                                           for x in md_line[modern_file.get("sample"):]]
                            write_to_file("\t".join([chrom,str(pos),ref,md_alt,archaic_geno]+md_geno)+"\n",filename)
                    else:
                        break
        return True
        #except:
           # return False
                        





    
        
## 先看看速度，考虑是否要将这些.split()提到前面去