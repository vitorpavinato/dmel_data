"""
This defines functions and classes to use in the annotation pipeline.
"""

import os
import subprocess
from pathlib import Path
from abc import ABC, abstractmethod
from dataclasses import dataclass
from simplify_snpeff import simplify_snpeff_default, simplify_snpeff_with_custom_annotation


@dataclass
class Annotation(ABC):
    """
    Abstract superclass for annotation classes.
    """
    input_file: str
    database: str
    output_folder: str

    def __post_init__(self):
        pass

    @abstractmethod
    def annotate(self):
        """
        Abstract method to perform the annotation.
        """
        pass


# SNPEff annotation subclass
@dataclass
class SNPEFFAnnotation(Annotation):
    """
    Class to annotate a vcf file with SNPEff.
    """
    SNPEFF = "/Users/tur92196/snpEff/snpEff.jar"
    SNPEFF_CONFIG = "/Users/tur92196/snpEff/snpEff.config"

    def __post_init__(self):
        # Define input folder and prefix file
        self.input_folder, filename = os.path.split(self.input_file)
        self.prefix, _ = filename.strip().split('.vcf')

        # Check if output folder exists; create it if it doesn't
        if not Path(self.output_folder).exists():
            Path(self.output_folder).mkdir(parents=True)

        self.output_file = f"{self.output_folder}/{self.prefix}_ann.vcf"
        self.csv_stats_file = f"{self.output_folder}/{self.prefix}_ann.csv"

    # Annotate a vcf using database
    def annotate(self):
        """
        Method to annotate the VCF file using SNPEff.
        """
        command = (
            f"java -jar {self.SNPEFF} ann -c {self.SNPEFF_CONFIG} -classic "
            f"-noStats -csvStats {self.csv_stats_file} -v {self.database} "
            f"{self.input_file} > {self.output_file}"
        )
        # Execute the command here
        subprocess.run(command, shell=True)

    # Annotate a vcf using only custom intervals
    def annotate_with_intervals(self, interval_file: str):
        """
        Method to annotate the VCF file using SNPEff
        using a custom interval file in BED format.
        """
        command = (
            f"java -jar {self.SNPEFF} ann -c {self.SNPEFF_CONFIG} -classic "
            f"-noStats -csvStats {self.csv_stats_file} -interval {interval_file} "
            f"-v {self.database} {self.input_file} > {self.output_file}"
        )
        # Execute the command here
        subprocess.run(command, shell=True)

    #@staticmethod
    def simplify_snpeff_annotations(self, custom_annotation: str = None):
        """
        Method to simplify the SNPEff annotations.
        """
        self.output_file_simplified = f"{self.output_folder}/{self.prefix}_ann_simplified.vcf"

        # Choose the simplification method based on the custom annotation
        if custom_annotation is not None:
            simplify_snpeff_with_custom_annotation(
                file=self.output_file, custom_annotation=custom_annotation
            )

        else:
            simplify_snpeff_default(file=self.output_file)


# SIFT4G annotation subclass
@dataclass
class SIFT4GAnnotation(Annotation):
    """
    Class to annotate a vcf file with SIFT4G.
    """

    SIFT4G = "/Users/tur92196/local/sift4g/SIFT4G_Annotator.jar"

    def __post_init__(self):
        # Check if output folder exists; create it if it doesn't
        if not Path(self.output_folder).exists():
            Path(self.output_folder).mkdir(parents=True)

    def annotate(self):
        """
        Method to annotate the VCF file using SIFT4G.
        """
        command = (
            f"java -jar {self.SIFT4G} -c -i {self.input_file} -d {self.database} -r {self.output_folder}"
        )
        # Execute the command here
        subprocess.run(command, shell=True)
