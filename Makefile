all:
	make -C samtools samtools
	cp -a samtools/samtools resources/usr/bin/samtools
	make -C samtools/bcftools bcftools
	cp -a samtools/bcftools/bcftools resources/usr/bin/bcftools
	cp -a samtools/bcftools/vcfutils.pl resources/usr/bin/vcfutils.pl

clean:
	rm -f resources/usr/bin/samtools resources/usr/bin/bcftools resources/usr/bin/vcfutils.pl
