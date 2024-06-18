import gzip


# 定义基类VCFFILE
class VCFFILE:
    # 传入参数：filename
    def __init__(self, filename: str):
        self._filename = filename
        self._maintext = self._get_maintext()

    # 判断file是否是gz压缩格式
    def _is_gzipped(self) -> bool:
        with open(self._filename, 'rb') as in_file:
            magic = in_file.read(2)
            return magic == b'\x1f\x8b'

    # 获取文本内容
    def _get_fulltext(self):
        if self._is_gzipped():
            with gzip.open(self._filename, 'rt') as in_file:
                full_text = in_file.readlines()
                return full_text
        else:
            with open(self._filename, 'rt') as in_file:
                full_text = in_file.readlines()
                return full_text

    # 去掉##注释
    def _get_maintext(self):
        # 使用header_count计数
        header_count = 0
        full_text = self._get_fulltext()
        for line in full_text:
            if line.startswith('##'):
                header_count += 1
            if line.startswith('#CHROM'):
                return full_text[header_count:]

# 定义现代人VCF子类


class MODERNVCF(VCFFILE):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._header_dict: dict = {}
        self._header_dict_setter()
        self._pos_pool: set = self._pos_pool_setter()

    # 定义获取文件内容的接口
    def get_text(self):
        return self._maintext

    # 获取样本header开始列，以FORMAT+1为基准,只用sample_index做header控制
    def _get_sampleindex(self) -> int:
        try:
            return self._maintext[0].strip().split("\t").index("FORMAT")+1
        except ValueError:
            print("Cannot find sample index，Check your header!\n")
            exit(1)
    # 将对应位置写入_header_dict

    def _header_dict_setter(self):
        header_line = self._maintext[0].strip().split("\t")
        try:
            self._header_dict.update({"chrom": header_line.index("#CHROM"),
                                      "pos": header_line.index("POS"),
                                      "ref": header_line.index("REF"),
                                      "alt": header_line.index("ALT"),
                                      "sample_index": self._get_sampleindex()})
        except ValueError:
            print("Header Error, Cannot find right index,Check your header in Modern VCF!\n")
            exit(1)
    # 读取对应属性的getter

    def _header_dict_getter(self, key: str):
        if key in self._header_dict.keys():
            return self._header_dict[key]
        else:
            print(f"This Property {key} is not present in modern VCF file.Please check!\n")
            exit(1)

    def get(self, key: str) -> int:
        return self._header_dict_getter(key)

    def _pos_pool_setter(self) -> set:
        _pos_pool = set()
        for line in self._maintext()[1:]:
            line = line.strip().split("\t")
            _pos_pool.add(line[1])
        return _pos_pool

    def _pos_pool_getter(self) -> set:
        return self._pos_pool

    def get_pospool(self) -> set:
        return self._pos_pool

    def get_samplecount(self) -> int:
        return len(self._maintext()[0].strip().split("\t")[self._get_sampleindex():])


class ARCHAICVCF(VCFFILE):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._header_dict = {}
        self._header_dict_setter()

    def get_text(self):
        return self._maintext()

    def _get_archaic_index(self) -> int:
        header_line = self._maintext()[0].strip().split("\t")
        try:
            return header_line.index("FORMAT")+1
        except ValueError:
            print("Cannot find archaic name index，Check your header!\n")
            exit(1)

    # 将对应位置写入_header_dict
    def _header_dict_setter(self):
        header_line = self._maintext()[0].strip().split("\t")
        try:
            self._header_dict.update({"chrom": header_line.index("#CHROM"),
                                      "pos": header_line.index("POS"),
                                      "ref": header_line.index("REF"),
                                      "alt": header_line.index("ALT"),
                                      "archaic_index": self._get_archaic_index(),
                                      "archaic_name": header_line[self._get_archaic_index()]})
        except ValueError:
            print("Header Error, Cannot find right index,Check your header in archaic VCF!\n")
            exit(1)
    # 读取对应属性的getter

    def _header_dict_getter(self, key: str):
        if key in self._header_dict.keys():
            return self._header_dict[key]
        else:
            print(f"This Property {key} is not present in archaic VCF file.Please check!\n")
            exit(1)

    def get(self, key: str) -> int:
        return self._header_dict_getter(key)
    # 获得染色体号

    def get_chrom(self):
        return self._maintext()[1].strip().split("\t")[self.get("chrom")]

    # 重写“+”号运算符，使用archaic+modern的方式实现merge
    def __add__(self, modern_file: MODERNVCF) -> bool:
        try:
            # 定义抽取指定列组合成新文件对象的函数:
            with gzip.open(f'{self.get("archaic_name")}_{self.get_chrom()}.gz', 'wt') as out_file:
                # 输出header
                out_file.write("\t".join(['chrom', 'pos', 'ref', 'alt', self.get("archaic_name")] +
                    modern_file.get_text()[0].strip().split("\t")[modern_file.get("sample_index"):]) + "\n")
                # 初始化现代人VCF文件指针
                seek_pointer = 1
                # 遍历古人中的每一行
                for line in self.get_text()[1:]:
                    line = line.strip().split("\t")
                    # 读取信息并存入，减少后续代码长度
                    chrom = line[self.get("chrom")]
                    pos = line[self.get("pos")]
                    ref = line[self.get("ref")]
                    alt = line[self.get("alt")]
                    # 0/0 以及./.的情况
                    try:
                        archaic_genotype = str(int(line[self.get("archaic_index")][0])
                                            + int(line[self.get("archaic_index")][2]))
                    except ValueError:
                        archaic_genotype = "0"
                    # 首先检查archaic_genotype非0且modern中没有出现的
                    if archaic_genotype != "0" and pos not in modern_file.get_pospool():
                        out_file.write("\t".join([chrom, pos, ref, alt, archaic_genotype]
                                                 + ["0"]*modern_file.get_samplecount()) + "\n")
                        continue
                    # 遍历现代人中的每一行，为减少消耗，在每一次大小的分界停下，取最后一个小位点作为下一次遍历开始的的文件指针
                    if seek_pointer < len(modern_file.get_text())-1:
                        for md_line in modern_file.get_text()[seek_pointer:]:
                            md_line = md_line.strip().split("\t")
                            md_pos = int(md_line[modern_file.get("pos")])
                            md_ref = md_line[modern_file.get("ref")]
                            md_alt = md_line[modern_file.get("alt")]
                            md_max_pos = int(modern_file.get_text()[-1].strip().split("\t")[modern_file.get("pos")])
                            if int(pos) > md_pos:
                                if int(pos) <= md_max_pos:
                                    seek_pointer += 1
                                break
                            elif int(pos) == md_pos:
                                if (len(md_ref) == 1 and len(md_alt) == 1 and md_ref == ref
                                        and (md_alt == alt or alt == '.')):
                                    md_genotype = [str(int(list(x)[0]) + int(list(x)[2]))
                                                   for x in md_line[modern_file.get("sample_index"):]]
                                    out_file.write("\t".join([chrom, pos, ref,
                                                              md_alt, archaic_genotype]+md_genotype)+"\n")
                            else:
                                break
            return True
        except:
            return False
