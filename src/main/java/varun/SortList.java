package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class SortList {

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

	public static <T extends RealType<T>> void listtoImage(PointSampleList<T> list,
			RandomAccessibleInterval<T> imageout) {

		final Cursor<T> testCursor = list.localizingCursor().copyCursor();
		final RandomAccess<T> imageCursor = imageout.randomAccess();

		while (testCursor.hasNext()) {

			testCursor.fwd();

			imageCursor.setPosition(testCursor);
			imageCursor.get().set(testCursor.get());

		}

	}

	public static <T extends RealType<T>> void split(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();
		if (list.dimension(direction) <= 1)
			return;

		else {
			if (direction == list.numDimensions())
				direction = 0;

			int otherdirection;

			otherdirection = direction + 1;
			if (direction == list.numDimensions() - 1)
				otherdirection = 0;
			if (direction >= 0 && direction < list.numDimensions() - 1)
				otherdirection = direction + 1;

			int meanIndex;

			final Cursor<T> listCursor = list.localizingCursor().copyCursor();

			PointSampleList<T> childA = new PointSampleList<T>(n);

			PointSampleList<T> childB = new PointSampleList<T>(n);

			meanIndex = (int) (list.min(direction) + ((list.max(direction) - list.min(direction)) / 2)
					+ list.dimension(direction) % 2);

			// In this loop I create the splitPlane at mean value in one
			// direction

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point splitPoint = new Point(n);

				long splitPlane = listCursor.getLongPosition(otherdirection);

				splitPoint.setPosition(meanIndex, direction);
				splitPoint.setPosition(splitPlane, otherdirection);

				Point cord = new Point(n);

				cord.setPosition(listCursor);
				if (listCursor.getLongPosition(direction) < splitPoint.getDoublePosition(direction)){

					childA.add(cord, listCursor.get().copy());
					
				}
				else{

					childB.add(cord, listCursor.get().copy());
				
				
				}
			}

			split(childA, direction);

			split(childB, direction);

			mergeList(list, childA, childB);
			 
		}

	}

	public static <T extends RealType<T>> void mergeList(PointSampleList<T> list, PointSampleList<T> listA,
			PointSampleList<T> listB) {

		
		int n= list.numDimensions();
		PointSampleList<T> tmplist = new PointSampleList<T>(n);
		
		
		

		
		Cursor<T> cursorA = listA.cursor();

		Cursor<T> cursorB = listB.cursor();
		
		Cursor<T> listcursor = tmplist.cursor();
		
		Cursor<T> copycursor = list.cursor();
		copycursor.fwd();
		while (copycursor.hasNext()) {
			copycursor.fwd();
			Point cord = new Point(n);

			cord.setPosition(copycursor);

			tmplist.add(cord, copycursor.get().copy());

		}

		
		
		listcursor.fwd();
		cursorA.fwd();
		cursorB.fwd();
		
		

		while (cursorA.hasNext() && cursorB.hasNext()) {

			
			Point cord = new Point(n);
			cord.setPosition(listcursor);

			if (cursorA.get().compareTo(cursorB.get()) < 0) {
				
				listcursor.get().set(cursorA.get().copy());
				listcursor.fwd();
				cursorA.fwd();
				tmplist.add(cord, listcursor.get().copy());
				

			} else {
				
				listcursor.get().set(cursorB.get().copy());
				listcursor.fwd();
				cursorB.fwd();
				
				tmplist.add(cord, listcursor.get().copy());

			}
			
		
			
		while (cursorA.hasNext() ) {
				
				
				listcursor.get().set(cursorA.get().copy());
				listcursor.fwd();
				cursorA.fwd();

				
				
	tmplist.add(cord, listcursor.get().copy());
				
				

			}
			
			while (cursorB.hasNext() ) {
				
				

				
				listcursor.get().set(cursorB.get().copy());
				listcursor.fwd();
				cursorB.fwd();
				tmplist.add(cord, listcursor.get().copy());
				
			}	
			
			

		}
		/*
		while(listcursor.hasNext()){
		
			listcursor.fwd();
			
			Point newpoint = new Point(n);
			newpoint.setPosition(listcursor);
			list.add(newpoint, listcursor.get());
		
		
		}*/
	
			
		//	System.out.println(list.size());
			
			
		

		

		

		}

		
		

	

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));
		Img<FloatType> imgout = new CellImgFactory<FloatType>().create(img, img.firstElement());
		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());
		int numPoints = 3;
		Random rnd = new Random(System.currentTimeMillis());
		FinalInterval interval = new FinalInterval(new long[] { 0, 20 });
		for (int i = 0; i < numPoints; ++i) {
			Point point = new Point(2);

			for (int d = 0; d < 2; ++d) {
				point.setPosition((int) (rnd.nextInt() * (interval.max(d) - interval.min(d)) + interval.min(d)), d);

				list.add(point, new FloatType(rnd.nextInt()));
			}
		}

		final ArrayList<FloatType> testList = new ArrayList<FloatType>();

		Cursor<FloatType> test = list.cursor();
		test.fwd();
		while (test.hasNext()) {
			test.fwd();
			testList.add(test.get().copy());

		}

		for (int i = 0; i < testList.size(); ++i) {

			System.out.println("This is initial list: " + testList.get(i));

		}

		// list = getList(img);

		split(list, 0);
		split(list, 1);

		final ArrayList<FloatType> testListtwo = new ArrayList<FloatType>();

		Cursor<FloatType> testtwo = list.cursor();
		testtwo.fwd();
		while (testtwo.hasNext()) {
			testtwo.fwd();
			testListtwo.add(testtwo.get().copy());

		}

		for (int i = 0; i < testListtwo.size(); ++i) {

			System.out.println("This should be sorted: " + testListtwo.get(i));

		}

		// listtoImage(list, imgout);

		// ImageJFunctions.show(imgout).setTitle("split along x and y");;

	}
}
