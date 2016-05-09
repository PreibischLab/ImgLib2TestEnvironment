package varun;

import java.io.File;

import com.sun.tools.javac.util.Pair;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.Watershed;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.LabelingType;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.integer.ShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.RealSum;
import net.imglib2.view.Views;
import util.ImgLib2Util;


public class Watershedding {

	
	@SuppressWarnings("deprecation")
	public static void OldWatersherImage(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> seedimg, int background) {

		int n = inputimg.numDimensions();
		long[] dimensions = new long[n];
		long[] positionseed = new long[n];

		for (int d = 0; d < n; ++d)
			dimensions[d] = inputimg.dimension(d);

		final NativeImgLabeling<Integer, ShortType> outputLabeling = new NativeImgLabeling<Integer, ShortType>(
				new ArrayImgFactory<ShortType>().create(inputimg, new ShortType()));

		final NativeImgLabeling<Integer, ShortType> seedLabeling = new NativeImgLabeling<Integer, ShortType>(
				new ArrayImgFactory<ShortType>().create(seedimg, new ShortType()));

		final Cursor<LabelingType<Integer>> seedlabelcursor = seedLabeling.localizingCursor();

		ImageJFunctions.show(inputimg).setTitle("Input image for watershedding");
		ImageJFunctions.show(seedimg).setTitle("Seed image for watershedding");

		final RandomAccess<FloatType> seedcursor = seedimg.randomAccess();

		// Fill the seedlabel image

		while (seedlabelcursor.hasNext()) {
			LabelingType<Integer> seedl = seedlabelcursor.next();
			seedlabelcursor.localize(positionseed);
			seedcursor.setPosition(positionseed);
			int seedLabel = (int) seedcursor.get().get();
			if (seedLabel == background)
				continue;
			seedl.setLabel(seedLabel);

		}

		final Watershed<FloatType, Integer> watershed = new Watershed<FloatType, Integer>();
		watershed.process();

		watershed.setSeeds(seedLabeling);
		watershed.setIntensityImage(inputimg);
		watershed.setStructuringElement(AllConnectedComponents.getStructuringElement(2));
		watershed.setOutputLabeling(outputLabeling);
		watershed.getResult();
		RandomAccessibleInterval<FloatType> outimg = new ArrayImgFactory<FloatType>().create(inputimg, new FloatType());
		Cursor<LabelingType<Integer>> labelCursor = outputLabeling.cursor();
		RandomAccess<FloatType> imageRA = outimg.randomAccess();

		// Go through the whole image again and again for every single label,
		// until no more label is found.
		int currentLabel = 0;
		boolean anythingFound = true;
		while (anythingFound) {
			anythingFound = false;
			labelCursor.reset();

			// Go through the whole image and add every pixel, that belongs to
			// the currently processed label
			int count = 0;
			while (labelCursor.hasNext()) {
				labelCursor.fwd();
				imageRA.setPosition(labelCursor);

				int i = (int) (imageRA.get().get());
				if (i == currentLabel) {
					anythingFound = true;
					count++;
				}
			}
			System.out.println("Number of input pixels in label " + currentLabel + ": " + count);
			currentLabel++;
			ImageJFunctions.show(outimg).setTitle("Watershed Images");
		}
		

	}
	

	public static enum InverseType {
		Straight, Inverse
	}
	// Thresholding a FlotType to convert to BitType
		public static void ThresholdingBit(RandomAccessibleInterval<FloatType> img,
				RandomAccessibleInterval<BitType> imgout, Float ThresholdValue) {

			final double[] backpos = new double[imgout.numDimensions()];
			final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

			final RandomAccess<BitType> outbound = imgout.randomAccess();

			while (bound.hasNext()) {

				bound.fwd();

				outbound.setPosition(bound);
	if (bound.get().get()>ThresholdValue){
//				if (bound.get().compareTo(ThresholdValue) > 0) {

					bound.localize(backpos);

					outbound.get().setReal(1);

				}

				else {

					outbound.get().setZero();

				}

			}
		}
		public static Pair<FloatType, FloatType> computeMinMaxIntensity(final RandomAccessibleInterval<FloatType> inputimg) {
			// create a cursor for the image (the order does not matter)
			final Cursor<FloatType> cursor = Views.iterable(inputimg).cursor();

			// initialize min and max with the first image value
			FloatType type = cursor.next();
			FloatType min = type.copy();
			FloatType max = type.copy();

			// loop over the rest of the data and determine min and max value
			while (cursor.hasNext()) {
				// we need this type more than once
				type = cursor.next();

				if (type.compareTo(min) < 0) {
					min.set(type);

				}

				if (type.compareTo(max) > 0) {
					max.set(type);

				}
			}
			Pair<FloatType, FloatType> pair = new Pair<FloatType, FloatType>(min, max);
			return pair;
		}

		// Automatic thresholding done on the Normalized input image
		// Algorithm: Get max and min intensity for an image and choose an initial
		// threshold value, T = (max-min)/2. This threshold value
		// segments image into two regions, get the mean pixel value for both the
		// regions (x1, x2)
		// then set the new threshold T_N = (x1 +x2)/2, segment initial image by
		// this value and repeat the process
		// till (T_N - T_{N+1}<epsilon) where epsilon is a small number say 1.0E-3
		public static Float AutomaticThresholding(RandomAccessibleInterval<FloatType> inputimg) {
			
			FloatType min = new FloatType();
			FloatType max = new FloatType();

			Float ThresholdNew, Thresholdupdate;

			Pair<FloatType, FloatType> pair = new Pair<FloatType, FloatType>(min, max);
			pair = computeMinMaxIntensity(inputimg);

			ThresholdNew = (pair.snd.get() - pair.fst.get()) / 2;

			// Get the new threshold value after segmenting the inputimage with thresholdnew
			Thresholdupdate = SegmentbyThresholding(Views.iterable(inputimg), ThresholdNew);

			while (true) {

				ThresholdNew = SegmentbyThresholding(Views.iterable(inputimg), Thresholdupdate);

				// Check if the new threshold value is close to the previous value
				if (Math.abs(Thresholdupdate - ThresholdNew) < 1.0E-2)
					break;
				Thresholdupdate = ThresholdNew;
			}
			

			return ThresholdNew;

		}

		// Segment image by thresholding, used to determine automatic thresholding
		// level
		public static Float SegmentbyThresholding(IterableInterval<FloatType> inputimg, Float Threshold) {

			int n = inputimg.numDimensions();
			Float ThresholdNew;
			PointSampleList<FloatType> listA = new PointSampleList<FloatType>(n);
			PointSampleList<FloatType> listB = new PointSampleList<FloatType>(n);
			Cursor<FloatType> cursor = inputimg.localizingCursor();
			while (cursor.hasNext()) {
				cursor.fwd();

				if (cursor.get().get() < Threshold) {
					Point newpointA = new Point(n);
					newpointA.setPosition(cursor);
					listA.add(newpointA, cursor.get().copy());
				} else {
					Point newpointB = new Point(n);
					newpointB.setPosition(cursor);
					listB.add(newpointB, cursor.get().copy());
				}
			}
			final RealSum realSumA = new RealSum();
			long countA = 0;

			for (final FloatType type : listA) {
				realSumA.add(type.getRealDouble());
				++countA;
			}

			final double sumA = realSumA.getSum() / countA;

			final RealSum realSumB = new RealSum();
			long countB = 0;

			for (final FloatType type : listB) {
				realSumB.add(type.getRealDouble());
				++countB;
			}

			final double sumB = realSumB.getSum() / countB;

			ThresholdNew = (float) (sumA + sumB) / 2;

			return ThresholdNew;

		}
	public static void DistanceTransformImage(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> outimg, final InverseType invtype) {
		int n = inputimg.numDimensions();

		final Img<BitType> bitimg = new ArrayImgFactory<BitType>().create(inputimg, new BitType());
		// make an empty list
		final RealPointSampleList<BitType> list = new RealPointSampleList<BitType>(n);

		final Float threshold = AutomaticThresholding(inputimg);
		ThresholdingBit(inputimg, bitimg, threshold);

		// cursor on the binary image
		final Cursor<BitType> cursor = bitimg.localizingCursor();

		// for every pixel that is 1, make a new RealPoint at that location
		while (cursor.hasNext())
			if (cursor.next().getInteger() == 1)
				list.add(new RealPoint(cursor), cursor.get());

		// build the KD-Tree from the list of points that == 1
		final KDTree<BitType> tree = new KDTree<BitType>(list);

		// Instantiate a nearest neighbor search on the tree (does not modifiy
		// the tree, just uses it)
		final NearestNeighborSearchOnKDTree<BitType> search = new NearestNeighborSearchOnKDTree<BitType>(tree);

		// randomaccess on the output
		final RandomAccess<FloatType> ranac = outimg.randomAccess();

		// reset cursor for the input (or make a new one)
		cursor.reset();

		// for every pixel of the binary image
		while (cursor.hasNext()) {
			cursor.fwd();

			// set the randomaccess to the same location
			ranac.setPosition(cursor);

			// if value == 0, look for the nearest 1-valued pixel
			if (cursor.get().getInteger() == 0) {
				// search the nearest 1 to the location of the cursor (the
				// current 0)
				search.search(cursor);

				// get the distance (the previous call could return that, this
				// for generality that it is two calls)
				switch (invtype) {

				case Straight:
					ranac.get().setReal(search.getDistance());
					break;
				case Inverse:
					ranac.get().setReal(-search.getDistance());
					break;

				}

			} else {
				// if value == 1, no need to search
				ranac.get().setZero();
			}
		}

	}
	
	// Finds and displays Local Maxima by constructing a 3*3*3.. local
		// neighbourhood
		public static RandomAccessibleInterval<FloatType> FindandDisplayLocalMaxima(RandomAccessibleInterval<FloatType> img,
				ImgFactory<FloatType> imageFactory) {

			// Create a new image for the output
			RandomAccessibleInterval<FloatType> output = imageFactory.create(img, new FloatType());

			// define an interval that is span number of pixel smaller on each side
			// in each dimension
			int span = 1;

			Interval interval = Intervals.expand(img, -span);

			// create a view on the source with this interval
			img = Views.interval(img, interval);

			// create a Cursor that iterates over the source and checks in a
			// 8-neighborhood
			// if it is a maxima
			final Cursor<FloatType> center = Views.iterable(img).cursor();

			// instantiate a RectangleShape to access rectangular local
			// neighborhoods

			final RectangleShape shape = new RectangleShape(span, true);

			// iterate over the set of neighborhoods in the image
			for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
				final FloatType centerValue = center.next();

				// keep this boolean true as long as no other value in the local
				// neighborhood
				// is smaller
				boolean isMaximum = true;

				// check if all pixels in the local neighborhood that are smaller
				for (final FloatType value : localNeighborhood) {
					// test if the center is smaller than the current pixel value
					if (centerValue.compareTo(value) < 0) {
						isMaximum = false;
						break;
					}
				}
				int n = img.numDimensions();
				double[] position = new double[n];
				if (isMaximum) {
					final RandomAccess<FloatType> outbound = output.randomAccess();
					outbound.setPosition(center);

					center.localize(position);
					
						outbound.get().set(1);
					
					

				}
			}

			return output;
		}
		public static void InvertInensityMap(RandomAccessibleInterval<FloatType> inputimg, FloatType minval, FloatType maxval){
	        // Normalize the input image
	 		Normalize.normalize(Views.iterable(inputimg), minval, maxval);
	 		// Now invert the normalization scale to get intensity inversion
			Normalize.normalize(Views.iterable(inputimg), maxval, minval);
		}
	public static void main(String[] args) {
	RandomAccessibleInterval<FloatType> biginputimg = ImgLib2Util
			.openAs32Bit(new File("src/main/resources/2015-01-14_Seeds-1.tiff"));
	
	
	new Normalize();
	FloatType minval = new FloatType(0);
	FloatType maxval = new FloatType(1);
	Normalize.normalize(Views.iterable(biginputimg), minval, maxval);
	final Img<FloatType> distimg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());

	DistanceTransformImage(biginputimg, distimg, InverseType.Inverse);

	ImageJFunctions.show(distimg).setTitle("DT inverted to extract seeds");

	RandomAccessibleInterval<FloatType> maximg = new ArrayImgFactory<FloatType>().create(biginputimg, new FloatType());

	final double[] sigma = { 0.5, 0.5 };

	maximg = FindandDisplayLocalMaxima(distimg, new ArrayImgFactory<FloatType>());

//	ImageJFunctions.show(maximg).setTitle("Seed Image for watershed");

	InvertInensityMap(distimg, minval, maxval);

//	ImageJFunctions.show(distimg).setTitle("DT image to perform watershed on");

	OldWatersherImage(distimg, maximg, 0);
	
	
	}
}
