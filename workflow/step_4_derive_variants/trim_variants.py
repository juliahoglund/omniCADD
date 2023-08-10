import sys, os, random
from optparse import OptionParser

# OptionParser for input. 
parser = OptionParser()
parser.add_option("-s", "--simu", dest="simu", help="Simulated variants", default=None)
parser.add_option("-d", "--derived", dest="derived", help="Derived variants", default=None)

(options, args) = parser.parse_args()

# Sorted list of files in directory
simu_file_list = []
derived_file_list = []
for fn in os.listdir(options.simu):
    if fn.startswith('snps_sim_variants_') and fn.endswith('_filtered.vcf'):
        simu_file_list.append(fn)
for fn in os.listdir(options.derived):
    if fn.startswith('derived_var_') and fn.endswith('_upper.vcf'):
        derived_file_list.append(fn)

simu_file_list = sorted(simu_file_list)
derived_file_list = sorted(derived_file_list)

print(simu_file_list)
print(derived_file_list)

# Count the total number of derived variants across all files
total_derived_variants = 0
for derived_fn in derived_file_list:
    with open(options.derived + derived_fn, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                total_derived_variants += 1
print(total_derived_variants)

# Shuffle the simulated variants and take as many as the total derived variants
for simu_fn in simu_file_list:
    output = open('trimmed_' + simu_fn, 'w')
    with open(options.simu + simu_fn, 'r') as f:
        all_lines = f.readlines()
        headers = [line for line in all_lines if line.startswith("#")]
        lines = [line for line in all_lines if not line.startswith("#")]
        if len(lines) < total_derived_variants:    
            print(f"Error: The simulated file {simu_fn} has fewer variants ({len(lines)}) than the total derived variants ({total_derived_variants}). Exiting.")
            sys.exit(1)

        # Print the first 5 lines before shuffling
        print("Before shuffling:")
        for l in lines[:5]:
            print(l.strip())

        random.shuffle(lines)

        # Print the first 5 lines after shuffling
        print("\nAfter shuffling:")
        for l in lines[:5]:
            print(l.strip())
        print("\n")
        for header in headers:
            print("header: "+ header)
            output.write(header)

        # Now, take the same number of lines as total_derived_variants
        for line in lines[:total_derived_variants]:
            output.write(line)

        output.close()

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_simulated_vcf_trimming.txt', 'x')
indication.close()
os.rename('./finished_simulated_vcf_trimming.txt', './output/finished_simulated_vcf_trimming.txt')
