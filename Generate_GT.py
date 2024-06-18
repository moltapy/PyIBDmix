from PramaParser import generate_gt_parser
from VCF_FILE_re import MODERNVCF
from VCF_FILE_re import ARCHAICVCF

# 获取解释器内的变量
gt_args = generate_gt_parser()

# 读入并创建实例
archaicvcf = ARCHAICVCF(gt_args.archaic)
modernvcf = MODERNVCF(gt_args.modern)


# 合并两个文件
merge = archaicvcf+modernvcf

# 条件判断

if merge:
    print(f"Generate GT Successfully For {gt_args.archaic} and {gt_args.modern}\n")
else:
    print(f"Error occurred In merge {gt_args.archaic} and {gt_args.modern}\n")
    exit(1)
