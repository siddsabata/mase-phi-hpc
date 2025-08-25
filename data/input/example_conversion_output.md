# SSM to DataFrame Conversion Example

## Input SSM Data (first 5 rows from data/ssm.txt):
```
id	gene	a	d	mu_r	mu_v
s0	GPR65_14_88477948_G>T	1287,1116,1262	2135,1233,1325	0.999	0.499
s1	GLIPR1L1_12_75728663_C>A	238,351,385	359,367,406	0.999	0.499
s2	DLG2_11_83544685_T>A	188,227,278	282,247,289	0.999	0.499
s3	LTA4H_12_96437132_A>C	360,332,332	525,362,341	0.999	0.499
s4	HYOU1_11_118925219_T>A	334,217,261	474,241,282	0.999	0.499
```

## Expected DataFrame Output (using samples 0 and 1):

| Hugo_Symbol | Reference_Allele | Allele | Chromosome | Start_Position | Variant_Frequencies_cf | Variant_Frequencies_st |
|-------------|------------------|--------|------------|----------------|------------------------|------------------------|
| GPR65       | G                | T      | 14         | 88477948       | 0.397                  | 0.095                  |
| GLIPR1L1    | C                | A      | 12         | 75728663       | 0.337                  | 0.044                  |
| DLG2        | T                | A      | 11         | 83544685       | 0.333                  | 0.081                  |
| LTA4H       | A                | C      | 12         | 96437132       | 0.314                  | 0.083                  |
| HYOU1       | T                | A      | 11         | 118925219      | 0.295                  | 0.100                  |

## VAF Calculations:

For each mutation, VAF = (total_depth - ref_count) / total_depth

**s0 (GPR65):**
- Sample 0: VAF = (2135 - 1287) / 2135 = 0.397
- Sample 1: VAF = (1233 - 1116) / 1233 = 0.095

**s1 (GLIPR1L1):**
- Sample 0: VAF = (359 - 238) / 359 = 0.337
- Sample 1: VAF = (367 - 351) / 367 = 0.044

## Key Features:

1. **Gene Parsing**: `GPR65_14_88477948_G>T` → Symbol: GPR65, Chr: 14, Pos: 88477948, Ref: G, Alt: T
2. **VAF Calculation**: Uses samples 0 and 1 as cf (circulating-free) and st (solid tissue) equivalents
3. **Filtering Preview**: All VAFs < 0.9, so no mutations would be filtered out
4. **Indexing**: Maintains original order (s0, s1, s2, ...) for compatibility with tree data

## Compatibility with Old Code:

This DataFrame format exactly matches what the old marker selection code expects:
- `Hugo_Symbol` → Gene symbol for naming
- `Reference_Allele`, `Allele` → For mutation annotation  
- `Chromosome`, `Start_Position` → For genomic coordinates in results
- `Variant_Frequencies_cf`, `Variant_Frequencies_st` → For VAF filtering (< 0.9)

The key advantage is that this conversion ensures the **same mutations in the same order** 
that would pass filtering, maintaining the critical indexing contract with the pre-computed 
tree distribution data. 