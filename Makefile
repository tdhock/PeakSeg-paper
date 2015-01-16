HOCKING-RIGAILL-chip-seq-paper.pdf: HOCKING-RIGAILL-chip-seq-paper.tex refs.bib figure-Segmentor-PeakSeg.png figure-dp-peaks-regression-dots.pdf figure-4samples-just-regions.png intervalRegression.pdf figure-interval-regression.tex figure-dp-peaks-train-2.png figure-regularized-all.pdf
	rm -f *.aux *.bbl
	pdflatex HOCKING-RIGAILL-chip-seq-paper
	bibtex HOCKING-RIGAILL-chip-seq-paper
	pdflatex HOCKING-RIGAILL-chip-seq-paper
	pdflatex HOCKING-RIGAILL-chip-seq-paper

figure-Segmentor-PeakSeg.png: figure-Segmentor-PeakSeg.R
	R --no-save < $<
figure-dp-peaks-regression-dots.pdf: figure-dp-peaks-regression-dots.R dp.peaks.regression.RData dp.peaks.baseline.RData regularized.all.RData
	R --no-save < $<
figure-dp-timings.pdf: figure-dp-timings.R dp.timings.RData
	R --no-save < $<

dp.peaks.baseline.RData: dp.peaks.baseline.R dp.peaks.sets.RData
	R --no-save < $<
dp.peaks.sets.RData: dp.peaks.sets.R dp.peaks.matrices.RData
	R --no-save < $<
dp.peaks.matrices.RData: dp.peaks.matrices.R dp.peaks.error.RData
	R --no-save < $<
dp.peaks.error.RData: dp.peaks.error.R dp.peaks.RData
	R --no-save < $<
dp.peaks.RData: dp.peaks.R dp.timings.RData
	R --no-save < $<
dp.timings.RData: dp.timings.R
	R --no-save < $<
dp.peaks.optimal.RData: dp.peaks.optimal.R dp.peaks.matrices.RData
	R --no-save < $<
dp.peaks.intervals.RData: dp.peaks.intervals.R dp.peaks.optimal.RData
	R --no-save < $<
dp.peaks.regression.RData: dp.peaks.regression.R dp.peaks.intervals.RData dp.peaks.sets.RData
	R --no-save < $<
dp.peaks.features.RData: dp.peaks.features.R dp.timings.RData
	R --no-save < $<

## For an interactive data viz comparing PeakSegDP against 2 baseline
## models on the benchmark data set.
figure-dp-peaks-interactive/index.html: figure-dp-peaks-interactive.R dp.peaks.interactive.RData dp.peaks.regression.RData
	R --no-save < $<
dp.peaks.interactive.RData: dp.peaks.interactive.R dp.peaks.regression.RData dp.peaks.baseline.RData
	R --no-save < $<

figure-interval-regression.tex: figure-interval-regression.R PeakSeg4samples.RData
	R --no-save < $<
intervalRegression.pdf: intervalRegression.tex figure-interval-regression.tex
	pdflatex intervalRegression
figure-4samples-just-regions.png: figure-PeakSeg-4samples.R PeakSeg4samples.RData
	R --no-save < $<
dp.peaks.train.RData: dp.peaks.train.R
	R --no-save < $<
PeakSeg4samples.RData: PeakSeg4samples.R dp.peaks.error.RData
	R --no-save < $<
figure-dp-peaks-train-2.png: figure-dp-peaks-train.R dp.peaks.train.RData
	R --no-save < $<

## regularized model for ICML paper.
regularized.RData: regularized.R dp.peaks.intervals.RData dp.peaks.sets.RData
	R --no-save < $<
figure-regularized.pdf: figure-regularized.R regularized.RData
	R --no-save < $<
regularized.one.RData: regularized.one.R dp.peaks.intervals.RData dp.peaks.sets.RData
	R --no-save < $<
regularized.all.RData: regularized.all.R dp.peaks.intervals.RData dp.peaks.sets.RData
	R --no-save < $<
figure-regularized-all.pdf: figure-regularized-all.R regularized.all.RData
	R --no-save < $<
