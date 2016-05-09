package varun;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.RealRandomAccess;
import net.imglib2.algorithm.gradient.PartialDerivative;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import varun.MyKDtree.Distance;
import varun.MyKDtree.Node;

public class KDtreejunk {
	
	/*
	 * 
	 	private static double gradientMagnitudeInDirection(double[] doublePos, double[] direction,
			RealRandomAccess<FloatType>[] gradientAccess) {
		double gradientmag = 0;
		for (int d = 0; d < gradientAccess.length; ++d) {
			gradientAccess[d].setPosition(doublePos);
			gradientmag += gradientAccess[d].get().getRealDouble() * direction[d];
		}
		return gradientmag;
	}
	 public static void CannyEdgeexp(RandomAccessibleInterval<FloatType> inputimg,
			RandomAccessibleInterval<FloatType> imgout, FloatType minval, FloatType maxval, double[] sigma,
			boolean Thresholding) {
		int n = inputimg.numDimensions();

		// n+1 th dimension is to index the partial derivatives of the n
		// dimensional image
		final long[] dim = new long[n + 1];
		for (int d = 0; d < n; ++d)
			dim[d] = inputimg.dimension(d);
		dim[n] = n;
		final Img<FloatType> gradients = new ArrayImgFactory<FloatType>().create(dim, inputimg.randomAccess().get());
		// RandomAccessibleInterval<FloatType> imgout = new
		// ArrayImgFactory<FloatType>().create(inputimg, new FloatType());

		final RandomAccess<FloatType> imgoutac = imgout.randomAccess();
		final long[] position = new long[n];

		final Float Minmagnitude = new Float(100);
		// GlobalThresholding.AutomaticThresholding(Views.iterable(gradients));
		System.out.println(Minmagnitude);
		// Compute partial derivatives of input in all dimension. This requires
		// a border of 1 pixel with respect to the input image
		final Interval gradientComputationInterval = Intervals.expand(inputimg, -1);
		for (int d = 0; d < n; ++d)
			PartialDerivative.gradientCentralDifference(inputimg,
					Views.interval(Views.hyperSlice(gradients, n, d), gradientComputationInterval), d);

		// To compute gradient maxima
		final Interval maximaComputationInterval = Intervals.expand(inputimg, -2);

		final long[] min = new long[n];
		final long[] max = new long[n];
		maximaComputationInterval.max(max);
		maximaComputationInterval.min(min);

		final long[] shiftback = new long[n];
		for (int d = 0; d < n; ++d)
			shiftback[d] = min[d] - max[d];

		final NLinearInterpolatorFactory<FloatType> interpolatorFactory = new NLinearInterpolatorFactory<FloatType>();
		@SuppressWarnings("unchecked")
		final RealRandomAccess<FloatType>[] gradientAccess = new RealRandomAccess[n];
		for (int d = 0; d < n; ++d)
			gradientAccess[d] = interpolatorFactory.create(Views.hyperSlice(gradients, n, d));

		final RandomAccess<FloatType> gradranac = gradients.randomAccess();

		for (int d = 0; d < n; ++d)
			gradranac.setPosition(min[d], d);
		gradranac.setPosition(0, n);

		final double direction[] = new double[n];
		final double doublePos[] = new double[n];

		final double minMagnitudeSquared = Minmagnitude * Minmagnitude;
		final long max0 = max[0];
		while (true) {
			// Get gradient direction and magnitude
			double gradientmag = 0;
			for (int d = 0; d < n; ++d) {
				final double partderiv = gradranac.get().getRealDouble();
				gradientmag += partderiv * partderiv;
				direction[d] = partderiv;
				gradranac.fwd(n);
			}
			gradranac.setPosition(0, n);
			if (gradientmag >= minMagnitudeSquared) {
				gradientmag = Math.sqrt(gradientmag);

				for (int d = 0; d < n; ++d) {
					direction[d] /= gradientmag;
					doublePos[d] = gradranac.getDoublePosition(d) + direction[d];
				}
				final double lighterMag = gradientMagnitudeInDirection(doublePos, direction, gradientAccess);
				if (gradientmag >= lighterMag) {
					for (int d = 0; d < n; ++d)
						doublePos[d] = gradranac.getDoublePosition(d) - direction[d];

					final double darkerMag = gradientMagnitudeInDirection(doublePos, direction, gradientAccess);

					if (gradientmag >= darkerMag) {
						// sub-pixel localization
						final double m = (darkerMag - lighterMag) / (2 * (darkerMag - 2 * gradientmag + lighterMag));
						for (int d = 0; d < n; ++d) {
							doublePos[d] = gradranac.getDoublePosition(d) + m * direction[d];
							position[d] = Math.round(doublePos[d]);
						}
						imgoutac.setPosition(position);
						imgoutac.get().setReal(gradientmag);
					}
				}
			}

			// move to next pixel
			if (gradranac.getLongPosition(0) == max0) {
				gradranac.move(shiftback[0], 0);
				// if ( n == 1 )
				// return edgels;
				for (int d = 1; d < n; ++d) {
					if (gradranac.getLongPosition(d) == max[d]) {
						gradranac.move(shiftback[d], d);
						// if ( d == n - 1 )
						// return edgels;
					} else {
						gradranac.fwd(d);
						break;
					}
				}
			} else
				gradranac.fwd(0);
		}

	}
	 
	// Find maxima only if the pixels lie along a given direction
	public static RandomAccessibleInterval<FloatType> FindDirectionalLocalMaxima(
			RandomAccessibleInterval<FloatType> img, ImgFactory<FloatType> imageFactory,
			final IntensityType setintensity, double[] sigma, double[] direction, Float val) {
		int n = direction.length;
		RandomAccessibleInterval<FloatType> output = imageFactory.create(img, new FloatType());
		// Construct a 5*5*5... local neighbourhood
		int span = 2;
		 double lambda=0;
		 double [] tmp = new double[direction.length];
		 
		Interval interval = Intervals.expand(img, -span);

		img = Views.interval(img, interval);
		
		double[] position = new double[n];
		final Cursor<FloatType> center = Views.iterable(img).cursor();

		final RectangleShape shape = new RectangleShape(span, true);

		for (final Neighborhood<FloatType> localNeighborhood : shape.neighborhoods(img)) {
			final FloatType centerValue = center.next();
			boolean isMaximum = true;
			final Cursor<FloatType> localcursor = localNeighborhood.localizingCursor();
		
			
			while (localcursor.hasNext()) {
				localcursor.fwd();
				for (int d = 0; d < n; ++d){
					if (direction[d]!=0){
						lambda = (center.getDoublePosition(d)-localcursor.getDoublePosition(d))/direction[d];
						break;
					}
					else{
					lambda = 0;
					}
				}
				
				double sum=0;
				for (int d = 0; d < n; ++d){
					if (direction[d]!=0){
					tmp[d] = (center.getDoublePosition(d)-localcursor.getDoublePosition(d))/direction[d];
					sum+=tmp[d];
					}
					else if(direction[d]==0){
						sum+=lambda;
					}
				}
				if (sum!=n*lambda && centerValue.compareTo(localcursor.get()) < 0 && centerValue.get() < val){
					System.out.println(lambda +" "+sum);
					isMaximum = false;
					break;
				}
				
			}
		
			
				if (isMaximum){
					
					final RandomAccess<FloatType> outbound = output.randomAccess();
				
					
						outbound.setPosition(center);
					center.localize(position);
					switch (setintensity) {

					case Original:
						outbound.get().set(center.get());
						break;

					case Gaussian:
						AddGaussian.addGaussian(output, position, sigma, false);
						break;

					case One:
						outbound.get().set(1);
						break;

					default:
						AddGaussian.addGaussian(output, position, sigma, false);
						break;

					}
					

				

			}

		}
		
		return output;
	}
	 	
	 	
	 	public void searchbyIndex(final RealLocalizable cursor, final int direction,
				int dirstartindex, int dirlastindex, int odirstartindex, int odirlastindex) throws FileNotFoundException {
			cursor.localize(Position);
			Bestdistsquared = Double.MAX_VALUE;
			Bestaxisdiffsquared = Double.MAX_VALUE;
			closestNodebyIndex(list, Xlist, Ylist, direction,
					 dirstartindex,  dirlastindex,  odirstartindex,  odirlastindex, Bestdistsquared);
		}
		
	 * 
	 public Node<T> makeNodebyIndex(PointSampleList<T> sortedlist, ArrayList<Point> Xlist, ArrayList<Point> Ylist,
				int direction, int dirstartindex, int dirlastindex, int odirstartindex, int odirlastindex) {

			final boolean directionchoice = direction == n - 1;

			final ArrayList<Point> XorYlist = directionchoice ? Ylist : Xlist;

			double[] point = new double[n];

			point = getMedianbyIndex(Xlist, Ylist, direction, dirstartindex, dirlastindex, odirstartindex,
					odirlastindex);

			final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
			final PointSampleList<T> RightTree = new PointSampleList<T>(n);
			final Cursor<T> listCursor = sortedlist.localizingCursor();
			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(n);

				cord.setPosition(listCursor);

				if (listCursor.getDoublePosition(direction) < point[direction])

					LeftTree.add(cord, listCursor.get());

				else

					RightTree.add(cord, listCursor.get());

			}

			Node<T> node = new Node<T>(point, direction, LeftTree, RightTree);

			return node;

		}

	private void closestNodebyIndex(final PointSampleList<T> list, final ArrayList<Point> Xlist,
				final ArrayList<Point> Ylist, final int direction,
				int dirstartindex, int dirlastindex, int odirstartindex, int odirlastindex, double olddist) throws FileNotFoundException {

			if (list.dimension(direction)<=2)   {

				return;
			}

			else  {

				final boolean directionchoice = direction == n - 1;
				final int otherdirection = directionchoice ? 0 : direction + 1;

				
				final Node<T> currentBest = makeNodebyIndex(list, Xlist, Ylist, direction,
						 dirstartindex, dirlastindex, odirstartindex, odirlastindex);

				double dist = 0;

				for (int d = 0; d < n; ++d) {

					dist += Math.pow((Position[d] - currentBest.nodePoint[d]), 2);
				}

				final double locationdiff = Position[currentBest.direction]
						- currentBest.nodePoint[currentBest.direction];

				final double axisdiff = locationdiff * locationdiff;

				final boolean leftbranchsearch = locationdiff < 0;

				if (dist <= Bestdistsquared) {

					Bestdistsquared = dist;
					finalnode = currentBest;

				}

				int newdirstartindex = directionchoice? odirstartindex:dirstartindex; // This index changes
				int newdirlastindex =  directionchoice? odirlastindex:dirlastindex; // This index changes
				int newodirstartindex = directionchoice? dirstartindex:odirstartindex; // This remains same
				int newodirlastindex =  directionchoice? dirlastindex:odirlastindex; // This remains same
				
				
				
			final	PointSampleList<T> searchBranch = leftbranchsearch ? currentBest.LeftTree : currentBest.RightTree;
				int newsearchdirstartindex = leftbranchsearch ? newdirstartindex : (newdirlastindex-newdirstartindex)/2;
				int newsearchdirlastindex = leftbranchsearch ? (newdirlastindex-newdirstartindex)/2:newdirlastindex; 

				final PointSampleList<T> nonsearchBranch = leftbranchsearch ? currentBest.RightTree
						: currentBest.LeftTree;
				int newnondirstartindex = leftbranchsearch ?  (newdirlastindex-newdirstartindex)/2:newdirstartindex ;
				int newnondirlastindex = leftbranchsearch ? newdirlastindex:(newdirlastindex-newdirstartindex)/2; 

			

				if(newsearchdirlastindex > newsearchdirstartindex && newnondirlastindex > newnondirstartindex  ){
					closestNodebyIndex(searchBranch, Xlist, Ylist, otherdirection,
							 newsearchdirstartindex, newsearchdirlastindex, newodirstartindex, newodirlastindex, dist);

		//		if (axisdiff <= Bestdistsquared)
			//		closestNodebyIndex(nonsearchBranch, Xlist, Ylist, otherdirection,
				//		 newnondirstartindex, newnondirlastindex, newodirstartindex, newodirlastindex, dist);
				}
			
			}

		}

		public double[] getMedianbyIndex(ArrayList<Point> Xlist, ArrayList<Point> Ylist, int direction,
				int dirstartindex, int dirlastindex, int odirstartindex, int odirlastindex) {

			final boolean directionchoice = direction == n - 1;
			final int otherdirection = directionchoice ? 0 : direction + 1;

			final ArrayList<Point> cordsort = directionchoice ? Ylist : Xlist;

			final ArrayList<Point> anticordsort = directionchoice ? Xlist : Ylist;

			final double[] medianPoint = new double[n];
			int medianindexA =  (dirlastindex - dirstartindex ) / 2;

			medianPoint[direction] = (cordsort.get(medianindexA).getDoublePosition(direction));

			int medianindexB =  (odirlastindex - odirstartindex ) / 2;

			medianPoint[otherdirection] = (anticordsort.get(medianindexB).getDoublePosition(otherdirection));

			return medianPoint;
		}

	/******
	 * Returns a root tree, I do this to initialize an ArrayList<Node<T>> in the
	 * main program which I overwrite later to include all the subtrees
	 ******/
/*
	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
	/*
		if (direction == list.numDimensions())
			direction = 0;
		if ((list.realMax(direction) - list.realMin(direction) + 1) <= 2)
			return null;

		else {

			double pivotElement;

			pivotElement = getMedian(list, direction);

			final PointSampleList<T> LeftTree = new PointSampleList<T>(n);
			final PointSampleList<T> RightTree = new PointSampleList<T>(n);

			final Cursor<T> listCursor = list.localizingCursor();

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point cord = new Point(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement)

					LeftTree.add(cord, listCursor.get().copy());

				else

					RightTree.add(cord, listCursor.get().copy());

			}

			Cursor<T> rightTreecursor = RightTree.localizingCursor();
			rightTreecursor.fwd();

			Node<T> node = new Node<T>(pivotElement, rightTreecursor, direction, LeftTree, RightTree);

			return node;

		}

	}

	

	/**********
	 * Starting the distance transform routine, done in a sort of non-Java
	 * fashion, by getting a needed object out of a void method by creating an
	 * empty ArrayList of the object and then overwritting it in the program.
	 * For a better tasting implementation see the overloaded version below
	 * which creates a NN object and uses getters and setters to get the correct
	 * distance transform.
	 * 
	 * @throws FileNotFoundException
	 **********/
/*
	public static <T extends RealType<T>> void distanceTransform(PointSampleList<BitType> list,
			PointSampleList<BitType> listzerosorones, RandomAccessibleInterval<T> imgout, final Distance dist)
					throws FileNotFoundException {

		Node<BitType> rootnode;

		rootnode = makeNode(list, 0);
		PrintStream out = new PrintStream(new FileOutputStream("KDtreemindist.txt"));
		System.setOut(out);

		double Bestdistsquared = Double.MAX_VALUE;

		final Cursor<BitType> zerooronelistcursor = listzerosorones.localizingCursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		double distance = 0;
		while (zerooronelistcursor.hasNext()) {
			zerooronelistcursor.fwd();
			double mindistance = Double.MAX_VALUE;
			double farmindistance = Double.MAX_VALUE;
			outbound.setPosition(zerooronelistcursor);

			ArrayList<Node<BitType>> nearnodelist = new ArrayList<Node<BitType>>(4);
			ArrayList<Node<BitType>> farnodelist = new ArrayList<Node<BitType>>(4);

			closestNode(zerooronelistcursor, rootnode, nearnodelist, Bestdistsquared);
			furtherNode(zerooronelistcursor, rootnode, farnodelist, Bestdistsquared);
			PointSampleList<BitType> singletree = combineTrees(nearnodelist.get(0));
			PointSampleList<BitType> singlefartree;
			if (farnodelist != null)
				singlefartree = combineTrees(farnodelist.get(0));
			else
				singlefartree = null;

			Cursor<BitType> singlecursor = singletree.cursor();

			while (singlecursor.hasNext()) {
				singlecursor.fwd();

				distance = dist.getDistance(zerooronelistcursor, singlecursor);
				mindistance = Math.min(distance, mindistance);

			}

			if (singlefartree != null) {

				Cursor<BitType> singlefarcursor = singlefartree.cursor();
				double fardistance = 0;
				while (singlefarcursor.hasNext()) {
					singlefarcursor.fwd();

					fardistance = dist.getDistance(zerooronelistcursor, singlefarcursor);
					farmindistance = Math.min(fardistance, farmindistance);

				}

			}
			final double actualmindistance = Math.min(mindistance, farmindistance);

			outbound.get().setReal(actualmindistance);

		}

		final Cursor<BitType> listcursor = list.localizingCursor();
		while (listcursor.hasNext()) {
			listcursor.fwd();
			outbound.setPosition(listcursor);
			outbound.get().setReal(0);

		}

	}
	/*
	/***********
	 * Returns the node closest to the testpoint stores as a single entry in the
	 * arrayList, done in a sort of "non-java" fashion
	 ***********/
/*
	public static <T extends RealType<T>> void closestNode(RealLocalizable testpoint, Node<T> Trees,
			ArrayList<Node<T>> list, double Bestdistsquared) {

		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - Trees.nodePoint[direction]);

		double dist = 0;
		for (int d = 0; d < n; ++d) {

			dist += (testpoint.getDoublePosition(d) - Trees.nodePoint[d])
					* (testpoint.getDoublePosition(d) - Trees.nodePoint[d]);
		}

		final Node<T> finalnode;

		final double mindistsquared = Math.min(dist, Bestdistsquared);

		if (dist < Bestdistsquared) {

			finalnode = Trees;
			Bestdistsquared = mindistsquared;

		}

		else
			finalnode = null;

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> searchBranch = leftbranchsearch ? Trees.LeftTree : Trees.RightTree;

		Node<T> nearnode;
		if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {

			nearnode = makeNode(searchBranch, otherdirection);
			if (nearnode != null)
				closestNode(testpoint, nearnode, list, mindistsquared);

		}

		if (finalnode != null)
			list.add(finalnode);

	}

	public static <T extends RealType<T>> void furtherNode(RealLocalizable testpoint, Node<T> Trees,
			ArrayList<Node<T>> list, double Bestdistsquared) {
		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else
			otherdirection = direction + 1;

		double locationdiff = (testpoint.getDoublePosition(direction) - Trees.nodePoint[direction]);

		double axisdiff = locationdiff * locationdiff;

		double dist = 0;
		for (int d = 0; d < n; ++d) {

			dist += (testpoint.getDoublePosition(d) - Trees.nodePoint[d])
					* (testpoint.getDoublePosition(d) - Trees.nodePoint[d]);
		}

		final Node<T> finalnode;

		final double mindistsquared = Math.min(dist, Bestdistsquared);

		if (dist < Bestdistsquared) {

			finalnode = Trees;

		}

		else
			finalnode = null;

		final boolean leftbranchsearch = locationdiff < 0;

		final PointSampleList<T> nonsearchBranch = leftbranchsearch ? Trees.RightTree : Trees.LeftTree;

		Node<T> farnode;

		if (axisdiff <= mindistsquared
				&& (nonsearchBranch.realMax(otherdirection) - nonsearchBranch.realMin(otherdirection) + 1) > 2) {
			farnode = makeNode(nonsearchBranch, otherdirection);
			if (farnode != null)
				furtherNode(testpoint, farnode, list, mindistsquared);
		}

		if (finalnode != null)
			list.add(finalnode);

	}
*/
}
