from setuptools import setup

setup(
    name="vatools",
    version="4.1.0",
    packages=["vatools"],
    entry_points={
        "console_scripts":[
            "vcf-expression-annotator = vatools.vcf_expression_annotator:main",
            "vcf-readcount-annotator = vatools.vcf_readcount_annotator:main",
            "vcf-info-annotator = vatools.vcf_info_annotator:main",
            "vcf-genotype-annotator = vatools.vcf_genotype_annotator:main",
            "vep-annotation-reporter = vatools.vep_annotation_reporter:main",
            "transform-split-values = vatools.transform_split_values:main",
        ]
    },
    install_requires=[
        'vcfpy',
        'pysam',
        'pandas',
        'gtfparse',
        'testfixtures',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5",

        'License :: OSI Approved :: MIT License',
    ],

    author = 'Susanna Kiwala, Chris Miller',
    author_email = "help@vatools.org",
    description = "A tool for annotating VCF files with expression and readcount data",
    license = "MIT License",
    keywords = "cancer sequencing variant variants gene expression readcounts VAF allele frequency FPKM TPM transcript VCF",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/vatools",
)
