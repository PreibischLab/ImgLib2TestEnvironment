package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import net.imglib2.PointSampleList;
import net.imglib2.RealPointSampleList;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Kdbranch {
	

	// To Create node from a List
	public static <T extends RealType<T>> Point createNode(PointSampleList<T> list, int direction) {

		// number of dimensions
		int n = list.cursor().numDimensions();

		int meanIndex;

		Point node = new Point(n);

		meanIndex = (int) list.dimension(direction) / 2;

		node.setPosition(meanIndex, direction);

		return node;

	}

	
	// Input a List and get two sublists split at the median
	public static <T extends RealType<T>> void getsubList(PointSampleList<T> list, int direction) {

		if (direction==list.numDimensions())
	     direction =0;
		int n = list.cursor().numDimensions();

		int meanIndex, endIndex;

		Point node = new Point(n);

		Point endpoint = new Point(n);

		meanIndex = (int) list.dimension(direction) / 2;
		endIndex = (int) list.dimension(direction);

		node.setPosition(meanIndex, direction);
		endpoint.setPosition(endIndex, direction);

		System.out.println("Initial Split point in this direction: " + node);
		System.out.println("Size of list in this direction: " + list.dimension(direction));

		PointSampleList<T> childA = new PointSampleList<T>(n);

		PointSampleList<T> childB = new PointSampleList<T>(n);

		final Cursor<T> listCursorA = list.localizingCursor();
		final Cursor<T> listCursorB = list.localizingCursor();

		while (listCursorA.hasNext()) {
			listCursorA.fwd();
			Point cord = new Point(n);

			cord.setPosition(listCursorA);
			if (listCursorA.getLongPosition(direction) < node.getDoublePosition(direction))
				childA.add(cord, listCursorA.get().copy());
		}
		
		if((meanIndex-1)>0 && childA.size()>2){	
			
			
			
			getsubList(childA, direction+1);
			
}


		while (listCursorB.hasNext()) {
			listCursorB.fwd();
			Point cordtwo = new Point(n);

			cordtwo.setPosition(listCursorB);
			if (listCursorB.getLongPosition(direction) >= node.getDoublePosition(direction))
				childB.add(cordtwo, listCursorB.get().copy());
		}
		
		if ((meanIndex + 1) < (list.size() - 1) && childB.size() > 2){		
		getsubList(childB, direction+1);
		
}
		

		/*
		 * for (int d = 0; d < n; ++d) {
		 * 
		 * meanIndex = (int) list.dimension(d) / 2; if ((meanIndex - 1) >= 0 &&
		 * childA.size() > 0) {
		 * 
		 * createNode(childA, direction + 1);
		 * 
		 * }
		 * 
		 * if ((meanIndex + 1) <= (list.size() - 1) && childB.size() > 0) {
		 * createNode(childB, direction + 1); }
		 * 
		 * }
		 */

	}

	public static <T extends RealType<T>> PointSampleList<T> getList(RandomAccessibleInterval<T> img) {

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
		PointSampleList<T> parent = new PointSampleList<T>(n);

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(n);

			cord.setPosition(first);

			parent.add(cord, first.get().copy());

		}

		return parent;

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		list = getList(img);

		getsubList(list, 0);

		// ImageJFunctions.show(img).setTitle("Original_Image");

		// ImageJFunctions.show(imgout).setTitle("Permuted_Image");

	}

}
