from setuptools import setup

setup(
    name="vcf-annotation-tools",
    version="0.0.4",
    packages=["vcf_annotation_tools"],
    entry_points={
        "console_scripts":[
            "vcf-expression-annotator = vcf_annotation_tools.vcf_expression_annotator:main",
            "vcf-readcount-annotator = vcf_annotation_tools.vcf_readcount_annotator:main",
        ]
    },
    install_requires=[
        'vcfpy',
        'pysam',
        'pandas',
        'gtfparse',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5",

        'License :: OSI Approved :: MIT License',
    ],

    author = "Susanna Kiwala",
    author_email = "ssiebert@wustl.edu",
    description = "A tool for annotating VCF files with expression and readcount data",
    license = "MIT License",
    keywords = "cancer sequencing variant variants gene expression readcounts VAF allele frequency FPKM TPM transcript VCF",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/vcf-annotation-tools",
)
