from setuptools import setup

setup(
    name="vcf-expression-encoder",
    version="0.0.1",
    entry_points={
        "console_scripts":[
            "vcf-expression-encoder = vcf_expression_encoder:main",
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

        "Programming Language :: Python :: 3.5"

        'License :: OSI Approved :: MIT License',
    ],

    author = "Susanna Kiwala",
    author_email = "ssiebert@wustl.edu",
    description = "A tool for annotating VCF files with expression data",
    license = "MIT License",
    keywords = "cancer sequencing variant variants gene expression transcript VCF",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/vcf-expression-encoder",
)
