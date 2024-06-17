import argparse


# Generate_GT软件的解释器
def generate_gt_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--archaic", help="Input the archaic VCF file,May contain multiple archaic \
    samples but only one will be utilized for IBDmix. Work for both compressed and uncompressed.", required=True)
    parser.add_argument("-m", "--modern", help="Input the modern VCF file. Work for both compressed and \
    uncompressed", required=True)
    args = parser.parse_args()
    return args


# IBDmix软件的解释器
def ibdmix_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", " --genotype", help="The input genotype file produces by\
     Generate_GT")
