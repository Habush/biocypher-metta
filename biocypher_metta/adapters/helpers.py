from inspect import getfullargspec
import hashlib
from math import log10, floor, isinf
from liftover import get_lifter

import hgvs.dataproviders.uta
from hgvs.easy import parser
from hgvs.extras.babelfish import Babelfish

ALLOWED_ASSEMBLIES = ['GRCh38']
_lifters = {}


def assembly_check(id_builder):
    def wrapper(*args, **kwargs):
        argspec = getfullargspec(id_builder)

        if 'assembly' in argspec.args:
            assembly_index = argspec.args.index('assembly')
            if assembly_index >= len(args):
                pass
            elif args[assembly_index] not in ALLOWED_ASSEMBLIES:
                raise ValueError('Assembly not supported')
        return id_builder(*args, *kwargs)

    return wrapper


@assembly_check
def build_variant_id(chr, pos_first_ref_base, ref_seq, alt_seq, assembly='GRCh38'):
    # pos_first_ref_base: 1-based position
    key = '{}_{}_{}_{}_{}'.format(str(chr).lower(), pos_first_ref_base, ref_seq, alt_seq, assembly)
    # return hashlib.sha256(key.encode()).hexdigest()
    return key

@assembly_check
def build_regulatory_region_id(chr, pos_start, pos_end, assembly='GRCh38'):
    # return '{}_{}_{}_{}_{}'.format(class_name, chr, pos_start, pos_end, assembly)
    return '{}_{}_{}_{}'.format(chr, pos_start, pos_end, assembly)


@assembly_check
def build_variant_id_from_hgvs(hgvs_id, validate=True, assembly='GRCh38'):
    # translate hgvs naming to vcf format e.g. NC_000003.12:g.183917980C>T -> 3_183917980_C_T
    if validate:  # use tools from hgvs, which corrects ref allele if it's wrong
        # got connection timed out error occasionally, could add a retry function
        hdp = hgvs.dataproviders.uta.connect()
        babelfish38 = Babelfish(hdp, assembly_name=assembly)
        try:
            chr, pos_start, ref, alt, type = babelfish38.hgvs_to_vcf(
                parser.parse(hgvs_id))
        except Exception as e:
            print(e)
            return None

        if type == 'sub' or type == 'delins':
            return build_variant_id(chr, pos_start + 1, ref[1:], alt[1:])
        else:
            return build_variant_id(chr, pos_start, ref, alt)

    # if no need to validate/query ref allele (e.g. single position substitutions) -> use regex match is quicker
    else:
        if hgvs_id.startswith('NC_'):
            chr = int(hgvs_id.split('.')[0].split('_')[1])
            if chr < 23:
                chr = str(chr)
            elif chr == 23:
                chr = 'X'
            elif chr == 24:
                chr = 'Y'
            else:
                print('Error: unsupported chromosome name.')
                return None

            pos_start = hgvs_id.split('.')[2].split('>')[0][:-1]
            if pos_start.isnumeric():
                ref = hgvs_id.split('.')[2].split('>')[0][-1]
                alt = hgvs_id.split('.')[2].split('>')[1]
                return build_variant_id(chr, pos_start, ref, alt)
            else:
                print('Error: wrong hgvs format.')
                return None
        else:
            print('Error: wrong hgvs format.')
            return None


# Arangodb converts a number to string if it can't be represented in signed 64-bit
# Using the approximation of a limit +/- 308 decimal points for 64 bits


def to_float(str):
    MAX_EXPONENT = 307

    number = float(str)

    if number == 0:
        return number

    if isinf(number) and number > 0:
        return float('1e307')

    if isinf(number) and number < 0:
        return float('1e-307')

    base10 = log10(abs(number))
    exponent = floor(base10)

    if abs(exponent) > MAX_EXPONENT:
        if exponent < 0:
            number = number * float(f'1e{abs(exponent) - MAX_EXPONENT}')
        else:
            number = number / float(f'1e{abs(exponent) - MAX_EXPONENT}')

    return number


def check_genomic_location(chr, start, end,
                           curr_chr, curr_start, curr_end):
    """
    Checks if the curr locations are within the specified locations (chr, start, end)
    If no chr is specified, then it returns True b/c that means we want to import all chromosomes
    Used when we want to filter the data imported from a file by location
    """
    if chr is None:  # import the data on all chromosomes
        return True
    else:  # filter by chromosome and (if specified) by location
        if chr != curr_chr:
            return False
        else:
            if start and end:
                if int(curr_start) >= start and int(curr_end) <= end:
                    return True
            elif start:
                if int(curr_start) >= start:
                    return True
            elif end:
                if int(curr_end) <= end:
                    return True
            else:
                return True
    return False


def convert_genome_reference(chr, pos, from_build='hg19', to_build='hg38'):
    """
    Convert a genomic coordinate from one reference build to another.

    Args:
        from_build (str): The reference build version to convert from (must be 'hg19' or 'hg38').
        to_build (str): The reference build version to convert to (must be 'hg19' or 'hg38', and different from `from_build`).
        chr (str): The chromosome identifier (e.g., 'chr1', 'chrX').
        pos (int): The genomic position on the chromosome.

    Returns:
        int: The converted genomic position in the target reference build, or None if the conversion fails.
    """
    if from_build not in ['hg19', 'hg38'] or to_build not in ['hg19', 'hg38'] or from_build == to_build:
        raise ValueError("Invalid reference build versions. 'from_build' and 'to_build' must be different and one of 'hg19' or 'hg38'.")

    lifter_key = f"{from_build}_{to_build}"

    # Initialize the lifter for the specified build conversion if not already cached
    if lifter_key not in _lifters:
        _lifters[lifter_key] = get_lifter(from_build, to_build)

    # Convert the chromosome identifier to a format compatible with the liftover library
    chr_no = chr.replace('chr', '').replace('ch', '')

    try:
        # Perform the liftover conversion using the cached lifter object
        converted = _lifters[lifter_key].query(chr_no, pos)[0][1]
        return int(converted)
    except:
        return None