package varun;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import varun.MyKDtree.Distance;
import varun.MyKDtree.Node;

public class KDtreejunk {
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

	/******* Returns the medianElement for input PointSampleList *******/

	public static <T extends RealType<T>> double getMedian(PointSampleList<T> list, int direction) {

		final Cursor<T> listCursor = list.localizingCursor();

		final ArrayList<Double> values = new ArrayList<Double>();

		while (listCursor.hasNext()) {
			listCursor.fwd();

			values.add(listCursor.getDoublePosition(direction));

		}

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
	
	private static <T> void nodetoList(final Node<T> node, final ArrayList<Node<T>> allnodes) {
		allnodes.add(node);
	}

	/******
	 * Returns a root tree, I do this to initialize an ArrayList<Node<T>> in the
	 * main program which I overwrite later to include all the subtrees
	 ******/

	public static <T extends RealType<T>> Node<T> makeNode(PointSampleList<T> list, int direction) {

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
	
	/********
	 * For a Node<T>, returns a single PointSampleList by combining the Left and
	 * Right Tree pairs into one
	 ***********/

	public static <T extends RealType<T>> PointSampleList<T> combineTrees(Node<T> Tree) {

		assert Tree.LeftTree.numDimensions() == Tree.RightTree.numDimensions();
		int n = Tree.LeftTree.numDimensions();

		PointSampleList<T> singleTree = new PointSampleList<T>(n);
		Cursor<T> treecursorA = Tree.LeftTree.cursor();
		Cursor<T> treecursorB = Tree.RightTree.cursor();

		while (treecursorA.hasNext()) {

			treecursorA.fwd();
			Point treepoint = new Point(n);

			treepoint.setPosition(treecursorA);

			singleTree.add(treepoint, treecursorA.get());

		}
		while (treecursorB.hasNext()) {

			treecursorB.fwd();
			Point treepoint = new Point(n);

			treepoint.setPosition(treecursorB);

			singleTree.add(treepoint, treecursorB.get().copy());

		}

		return singleTree;

	}
	
	/***********
	 * Returns the node closest to the testpoint stores as a single entry in the
	 * arrayList, done in a sort of "non-java" fashion
	 ***********/

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

}
