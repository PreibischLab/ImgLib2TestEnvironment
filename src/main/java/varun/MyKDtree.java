package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import net.imglib2.RealPointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MyKDtree {

	/********* For an input image returns a RealPointSampleList ********/
	public static <T extends RealType<T>> RealPointSampleList<T> getList(RandomAccessibleInterval<T> img) {

		final RandomAccessible<T> infinite = Views.extendZero(img);

		final int n = img.numDimensions();
		long min[] = new long[n];
		long max[] = new long[n];

		for (int d = 0; d < n; ++d) {

			min[d] = img.min(d);
			max[d] = img.max(d);

		}

		FinalInterval interval = new FinalInterval(min, max);

		final IterableInterval<T> imgav = Views.interval(infinite, interval);

		final Cursor<T> first = imgav.cursor();

		// A point sample list with coordinates declared and initialized.
		RealPointSampleList<T> parent = new RealPointSampleList<T>(n);

		while (first.hasNext()) {
			first.fwd();
			RealPoint cord = new RealPoint(n);

			cord.setPosition(first);

			parent.add(cord, first.get().copy());

		}

		return parent;

	}

	/*********
	 * Starting the methods which sort an Arraylist of co-ordinates using
	 * Merge-Sort algorithm
	 *********/

	public static <T extends RealType<T>> void split(ArrayList<Double> coordinateList, int direction) {

		if (coordinateList.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) coordinateList.size() / 2;

			Iterator<Double> iterator = coordinateList.iterator();

			final ArrayList<Double> childA = new ArrayList<Double>((int) coordinateList.size() / 2);

			final ArrayList<Double> childB = new ArrayList<Double>(
					(int) (coordinateList.size() / 2 + coordinateList.size() % 2));

			int xindex = 0;

			while (iterator.hasNext()) {
				iterator.next();

				if (xindex < splitIndex)
					childA.add(coordinateList.get(xindex));

				else

					childB.add(coordinateList.get(xindex));

				xindex++;

			}

			split(childA, direction);

			split(childB, direction);

			mergeListValue(coordinateList, childA, childB);

		}

	}

	/// ***** Returns a sorted list *********////
	public static <T extends RealType<T>> void mergeListValue(ArrayList<Double> sortedlist, ArrayList<Double> listA,
			ArrayList<Double> listB) {

		int i = 0, j = 0, k = 0;

		while (i < listA.size() && j < listB.size()) {

			if (listA.get(i) < listB.get(j)) {

				sortedlist.set(k, listA.get(i));

				++i;
				++k;
			}

			else {

				sortedlist.set(k, listB.get(j));

				++j;
				++k;

			}

		}

		while (i < listA.size()) {
			sortedlist.set(k, listA.get(i));
			++i;
			++k;

		}

		while (j < listB.size()) {
			sortedlist.set(k, listB.get(j));
			++j;
			++k;

		}

	}
	/******** End of the Merge-Sort routine for Arraylist *********/

	/******* Returns the medianElement for input RealPointSampleList *******/

	public static <T extends RealType<T>> double getMedian(RealPointSampleList<T> list, int direction) {

		final RealCursor<T> listCursor = list.localizingCursor();

		final ArrayList<Double> values = new ArrayList<Double>();

		while (listCursor.hasNext()) {
			listCursor.fwd();

			values.add(listCursor.getDoublePosition(direction));

		}

		// Collections.sort(values);

		split(values, direction); // Since the list is sorted I only have to get
									// the value at the middle of the list

		int startindex = 0;
		int lastindex = values.size() - 1;

		// Size of a list is lastindex-startindex+1.

		int[] medianindex = new int[2];
		int medianIndexA, medianIndexB;

		if ((lastindex - startindex + 1) % 2 == 1) {
			medianIndexA = startindex + (lastindex - startindex + 1 - (lastindex - startindex + 1) % 2) / 2;
			medianIndexB = medianIndexA;
		}

		else {

			medianIndexA = startindex + (lastindex - startindex + 1) / 2 - 1;
			medianIndexB = medianIndexA + 1;
		}

		medianindex[0] = medianIndexA;

		medianindex[1] = medianIndexB;

		double medianElement = 0.0;

		medianElement = 0.5 * (values.get(medianindex[0]) + values.get(medianindex[1]));

		return medianElement;

	}

	/********
	 * Constructor for the object Node that contains the Value at which a list
	 * is split up, the two split lists and the direction of the split
	 *********/

	private static class searchNode<T> {

		private final int n;

		private final double medianValue;

		private final int direction;

		private final RealPointSampleList<T> searchBranch;

		public searchNode(final double medianValue, final int direction, final RealPointSampleList<T> searchBranch) {

			this.n = searchBranch.numDimensions();
			this.medianValue = medianValue;
			this.direction = direction;
			this.searchBranch = searchBranch;

		}

		public int getDirection() {
			return direction;
		}

		public double getMedianValue() {
			return medianValue;
		}

		public RealPointSampleList<T> getSearchBranch() {
			return searchBranch;
		}

		public int getNumdimensions() {
			return n;
		}

	}

	private static class Node<T> {

		public final int n;

		public final double medianValue;

		public final int direction;

		public final RealPointSampleList<T> LeftTree;

		public final RealPointSampleList<T> RightTree;

		public Node(final double medianValue, final int direction, final RealPointSampleList<T> LeftTree,
				final RealPointSampleList<T> RightTree) {
			assert LeftTree.numDimensions() == RightTree.numDimensions();
			this.n = LeftTree.numDimensions();
			this.medianValue = medianValue;
			this.direction = direction;
			this.LeftTree = LeftTree;
			this.RightTree = RightTree;
		}

		public int getnumDimensions() {
			return n;
		}

		public double getMedianValue() {
			return medianValue;

		}

		public int getDirection() {
			return direction;
		}

		public RealPointSampleList<T> getLeftTree() {
			return LeftTree;
		}

		public RealPointSampleList<T> getRightTree() {
			return RightTree;
		}

	}

	private static <T> void nodetoList(final Node<T> node, final ArrayList<Node<T>> allnodes) {
		allnodes.add(node);
	}

	private static <T> void searchNodetoList(final searchNode<T> searchnode, final ArrayList<searchNode<T>> allnodes) {
		allnodes.add(searchnode);
	}

	/******
	 * Returns a root tree, I do this to initialize an ArrayList<Node<T>> in the
	 * main program which I overwrite later to include all the subtrees (Clever
	 * or Dangerous?)
	 ******/

	public static <T extends RealType<T>> Node<T> makeNode(RealPointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == list.numDimensions())
			direction = 0;
		if ((list.realMax(direction) - list.realMin(direction) + 1) <= 2)
			return null;

		else {

			double pivotElement;

			pivotElement = getMedian(list, direction);

			final RealPointSampleList<T> LeftTree = new RealPointSampleList<T>(n);
			final RealPointSampleList<T> RightTree = new RealPointSampleList<T>(n);

			final RealCursor<T> listCursor = list.localizingCursor();

			while (listCursor.hasNext()) {

				listCursor.fwd();

				RealPoint cord = new RealPoint(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement) {

					LeftTree.add(cord, listCursor.get().copy());

				} else {

					RightTree.add(cord, listCursor.get().copy());

				}

			}

			Node<T> node = new Node<T>(pivotElement, direction, LeftTree, RightTree);

			return node;

		}

	}

	/*********
	 * Here I return an Arraylist of Node<T> type by a clever trick, I give in
	 * an ArrayList of Node<T> containing only the rootNode and in the course of
	 * the routine below overwrite that list with the nodes of all the subtrees.
	 * The list then contains the object Node<T> for all the subtrees including
	 * the rootTree. Is this clever or dangerous way to do it?
	 ********/
	public static <T extends RealType<T>> void getTree(RealPointSampleList<T> list, ArrayList<Node<T>> allnodes,
			int direction) {

		int n = list.numDimensions();
		/****
		 * To ward against running over the dimensionality, creating some local
		 * restrictions on the global variable direction
		 ****/
		if (direction == n)
			direction = 0;
		if ((list.realMax(direction) - list.realMin(direction) + 1) <= 2)
			return;

		else {

			// ArrayList<Node<T>> allnodes = new ArrayList<Node<T>>();

			double pivotElement;

			pivotElement = getMedian(list, direction);

			final RealPointSampleList<T> LeftTree = new RealPointSampleList<T>(n);
			final RealPointSampleList<T> RightTree = new RealPointSampleList<T>(n);

			final RealCursor<T> listCursor = list.localizingCursor();

			while (listCursor.hasNext()) {

				listCursor.fwd();

				RealPoint cord = new RealPoint(listCursor);

				if (listCursor.getDoublePosition(direction) < pivotElement) {

					LeftTree.add(cord, listCursor.get().copy());

				} else if (listCursor.getDoublePosition(direction) >= pivotElement) {

					RightTree.add(cord, listCursor.get().copy());

				}

			}

			Node<T> node = new Node<T>(pivotElement, direction, LeftTree, RightTree);

			int otherdirection = direction + 1;

			if (otherdirection == n)
				otherdirection = 0;

			getTree(LeftTree, allnodes, otherdirection);

			getTree(RightTree, allnodes, otherdirection);

			nodetoList(node, allnodes);

		}

	}

	/***********
	 * Returns the complete search path (splitNodes, direction of split and the
	 * search branches) giving all the nearest neighbours of a point, also
	 * stored are the farther neighbours of the testpoint,
	 * 
	 * index 0 of the nodeList stores the closest node and the search branch and
	 * lastindex stores the left or right side or the Root Tree where the first
	 * split happened to search for the point.
	 * 
	 * index 0 of the farnodeList stores the least farthest node to the search
	 * point and the lastindex stores the other side of the RootTree which
	 * should really be far far away from the given point.
	 ***********/

	public static <T extends RealType<T>> void closestNode(double[] testpoint, Node<T> Trees,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList) {

		int direction = Trees.direction;

		int n = Trees.getnumDimensions();
		int otherdirection;

		if (direction == n - 1)

			otherdirection = 0;

		else

			otherdirection = direction + 1;

		double locationdiff = (testpoint[direction] - Trees.getMedianValue());

		final boolean leftbranchsearch = locationdiff < 0;

		final RealPointSampleList<T> searchBranch = leftbranchsearch ? Trees.LeftTree : Trees.RightTree;
		final RealPointSampleList<T> nonsearchBranch = leftbranchsearch ? Trees.RightTree : Trees.LeftTree;

		Node<T> nextnode;
		searchNode<T> searchnode;
		searchNode<T> nonsearchnode;

		if ((searchBranch.realMax(otherdirection) - searchBranch.realMin(otherdirection) + 1) > 2) {
			nextnode = makeNode(searchBranch, otherdirection);
			searchnode = new searchNode<T>(nextnode.medianValue, otherdirection, searchBranch);
			nonsearchnode = new searchNode<T>(nextnode.medianValue, otherdirection, nonsearchBranch);
		}

		else {

			nextnode = Trees;
			searchnode = new searchNode<T>(nextnode.medianValue, otherdirection, searchBranch);
			nonsearchnode = new searchNode<T>(nextnode.medianValue, otherdirection, nonsearchBranch);

		}

		if (nextnode != Trees) {

			closestNode(testpoint, nextnode, nodeList, farnodeList);
		}

		searchNodetoList(searchnode, nodeList);
		searchNodetoList(nonsearchnode, farnodeList);

	}

	public static <T extends RealType<T>> double volumeHypercube(RealPointSampleList<T> list, int dimensions) {
		double vol = list.realMax(dimensions) - list.realMin(dimensions) + 1;
		// dimensions (of space) = list.numDimensions() - 1;

		for (int d = dimensions - 1; d >= 0; --d)
			vol *= list.realMax(d) - list.realMin(d) + 1;

		return vol;

	}

	public static <T extends RealType<T>> double NearestNeighbourSearch(double[] testpoint,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList, RealPointSampleList<T> list) {

		int n = list.numDimensions();

		int dimensions = n - 1;

		double Volume = volumeHypercube(list, dimensions);

		int depth = nodeList.size();

		double constantfactor = Math.sqrt(dimensions) * Math.pow(Volume, 1.0 / dimensions);

		double smallRadius = constantfactor / (Math.pow(2, depth / dimensions + 1));

		double bigRadius = constantfactor / (Math.pow(2, (depth - 1) / dimensions + 1));

		System.out.println(smallRadius);
		System.out.println(bigRadius);

		return smallRadius;

	}

	public static <T extends RealType<T>> Pair<Double , searchNode<T> > NearestNeighbourSearch(double[] testpoint,
			ArrayList<searchNode<T>> nodeList, ArrayList<searchNode<T>> farnodeList, final Distance dist) {

		double mindistance;
		double bestdistance = Double.MAX_VALUE;
		double secondbestdistance = Double.MAX_VALUE;
		double thirdbestdistance = Double.MAX_VALUE;

		searchNode<T> finalNode, semifinalNode;

		final RealCursor<T> listcursor = nodeList.get(0).getSearchBranch().localizingCursor();

		finalNode = nodeList.get(0);

		while (listcursor.hasNext()) {
			listcursor.fwd();
			mindistance = dist.getDistance(listcursor, testpoint);

			bestdistance = Math.min(mindistance, bestdistance);

		}
		for (int index = 1; index < nodeList.size() - 1; ++index) {

			final RealCursor<T> cursor = nodeList.get(index).getSearchBranch().localizingCursor();

			while (cursor.hasNext()) {
				cursor.fwd();
				mindistance = dist.getDistance(cursor, testpoint);
				secondbestdistance = Math.min(mindistance, secondbestdistance);

				if (secondbestdistance < bestdistance) {
					bestdistance = secondbestdistance;
					finalNode = nodeList.get(index);
					continue;
				}

				else {

					break;
				}
			}
		}

		final RealCursor<T> farlistcursor = farnodeList.get(0).getSearchBranch().localizingCursor();

		semifinalNode = farnodeList.get(0);

		while (farlistcursor.hasNext()) {
			farlistcursor.fwd();
			mindistance = dist.getDistance(farlistcursor, testpoint);

			thirdbestdistance = Math.min(mindistance, bestdistance);

			if (thirdbestdistance < bestdistance) {
				bestdistance = thirdbestdistance;
				finalNode = semifinalNode;
continue;
			}
		
			else
				continue;
			
			
		}

		

		for (int index = 1; index < farnodeList.size() - 1; ++index) {

			final RealCursor<T> farcursor = farnodeList.get(index).getSearchBranch().localizingCursor();

			while (farcursor.hasNext()) {
				farcursor.fwd();
				mindistance = dist.getDistance(farcursor, testpoint);
				thirdbestdistance = Math.min(mindistance, thirdbestdistance);

				if (thirdbestdistance < bestdistance) {
					bestdistance = thirdbestdistance;
					finalNode = farnodeList.get(index);

					continue;
				}

				else {

					break;
				}
			}
		}

		Pair<Double , searchNode<T> > pair = new ValuePair<Double, searchNode<T>>(bestdistance, finalNode);
		
		return pair;

	}

	

	/**********
	 * Starting the distance transform routine 
	 **********/

	
	public static  void distanceTransform(RealPointSampleList<BitType> list){
		
		int n = list.numDimensions();
		
		RealPointSampleList<BitType> bitlist;
		Node<BitType> rootnode;
		ArrayList<searchNode<BitType>> searchnodes = new ArrayList<searchNode<BitType>>();

		ArrayList<searchNode<BitType>> nonsearchnodes = new ArrayList<searchNode<BitType>>();
		
		rootnode = makeNode(list, 0);
		
		ArrayList<Node<BitType>> allnodes = new ArrayList<Node<BitType>>();
		
		
		
		RealCursor<BitType> listcursor = list.localizingCursor();
		
		while(listcursor.hasNext()){
			listcursor.fwd();
			
			if (listcursor.get().getInteger() == 0){
				
				// Start finding the nearest neighbours for all the points having 0 value
				
				double[] testpoint= new double[n] ;
				
				for (int d=0; d < n; ++d)
					testpoint[d] = listcursor.getDoublePosition(d);
				
				
				closestNode(testpoint, rootnode, searchnodes, nonsearchnodes);
				
				
				
				
				
			}
			
			
			
		}
		
		
		
		
		
		
	}
	
	
	
	public interface Distance {

		double getDistance(RealLocalizable cursor1, RealLocalizable cursor2);

		<T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint);

	}

	public static class EucledianDistance implements Distance {
		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {

			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.pow((cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d)), 2);

			}

			return Math.sqrt(distance);
		}

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += (testpoint[d] - listcursor.getDoublePosition(d))
						* (testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	public static class MannhattanDistance implements Distance {

		public double getDistance(RealLocalizable cursor1, RealLocalizable cursor2) {
			double distance = 0.0;

			for (int d = 0; d < cursor2.numDimensions(); ++d) {

				distance += Math.abs(cursor1.getDoublePosition(d) - cursor2.getDoublePosition(d));

			}

			return distance;
		}

		public <T extends RealType<T>> double getDistance(RealLocalizable listcursor, double[] testpoint) {
			double distance = 0.0;
			int n = listcursor.numDimensions();

			for (int d = 0; d < n; ++d)
				distance += Math.abs(testpoint[d] - listcursor.getDoublePosition(d));

			return Math.sqrt(distance);
		}

	}

	/************
	 * Creating a bitType image from an image of type T by doing thresholding
	 * (not used currently)
	 ************/

	public static <T extends RealType<T>> void createBitimage(Img<T> img, Img<BitType> imgout, T ThresholdValue) {

		final Cursor<T> bound = img.localizingCursor();

		final RandomAccess<BitType> outbound = imgout.randomAccess();

		while (bound.hasNext()) {

			bound.fwd();

			outbound.setPosition(bound);

			if (bound.get().compareTo(ThresholdValue) > 0) {

				outbound.get().setOne();

			}

			else {

				outbound.get().setZero();

			}

		}
	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		final Img<BitType> imgout = new ArrayImgFactory<BitType>().create(img, new BitType());
		
		FloatType val = new FloatType(200);
		
		createBitimage(img, imgout, val);
		
		RealPointSampleList<FloatType> list = new RealPointSampleList<FloatType>(img.numDimensions());
		RealPointSampleList<BitType> bitlist = new RealPointSampleList<BitType>(imgout.numDimensions());
		
		// Make a list by setting an appropriate
		// interval on the image.

		RandomAccessibleInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 5, 5 });

		list = getList(img);
		bitlist= getList(imgout);

		
		final RealCursor<FloatType> first = list.cursor();

		while (first.hasNext()) {
			first.fwd();

			// System.out.println(" list X cor:" + first.getDoublePosition(0));
			// System.out.println(" list Y cor:" + first.getDoublePosition(1));
			// System.out.println(" list: " + first.get());

		}

		int n = list.numDimensions();

		/********** Starting the KD-Tree creation *********/

		Node<FloatType> rootnode, finalnode;

		rootnode = makeNode(list, 0);

		ArrayList<Node<FloatType>> allnodes = new ArrayList<Node<FloatType>>();

		ArrayList<searchNode<FloatType>> searchnodes = new ArrayList<searchNode<FloatType>>();

		ArrayList<searchNode<FloatType>> nonsearchnodes = new ArrayList<searchNode<FloatType>>();

		// Make a KD-tree by splitting along the X direction first and then the
		// Y direction and so on until there is no more splitting possible
		getTree(list, allnodes, 0);

		
		/******** Make a test point and search for the closest node *********/

		double[] testpoint = new double[2];

		testpoint[0] = 0.4;
		testpoint[1] = 1.2;

		closestNode(testpoint, rootnode, searchnodes, nonsearchnodes);

		
		
		Pair< Double , searchNode<FloatType> > pair;
		
		pair = NearestNeighbourSearch(testpoint, searchnodes, nonsearchnodes, new EucledianDistance());

		System.out.println("Distance to Nearest Neighbour: " +pair.getA());
		
		
		final RealCursor<FloatType> testcursor = pair.getB().searchBranch.localizingCursor();
while(testcursor.hasNext()){
	
	testcursor.fwd();
	
	System.out.println("X-cord of points on nearest branch: " +testcursor.getDoublePosition(0));
	System.out.println("Y-cord of points on nearest branch: " +testcursor.getDoublePosition(1));
	System.out.println("PxVal of points on nearest branch:  "+testcursor.get());
	
}
		
	}

}
