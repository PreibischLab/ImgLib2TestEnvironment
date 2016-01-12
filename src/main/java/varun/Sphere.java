package varun;

import java.util.Random;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.Sampler;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class Sphere {

	public static <T extends RealType<T>> void drawMulti(final RandomAccessibleInterval<T> img, Point newcenter,
			long newradius) {

		HyperSphere<T> bigSphere = new HyperSphere<T>(img, newcenter, newradius);

		Random ram = new Random(10);

		Cursor<T> bigc = bigSphere.cursor();

		while (bigc.hasNext()) {

			bigc.fwd();

			if (ram.nextDouble() < 0.01) {

				double smvar = ram.nextDouble();
				long newradiustwo = Math.round(ram.nextDouble() * 10);

				HyperSphere<T> smallSphere = new HyperSphere<T>(img, bigc, newradiustwo);
				Cursor<T> smallc = smallSphere.cursor();

				while (smallc.hasNext()) {

					smallc.fwd();
					if (smvar > smallc.get().getRealDouble()) {

						T value = smallc.get(); //Gives the value at that pixel
						
						value.setReal(smvar); // Sets the value

					}

				}
			}

		}

	}

	public static <T extends RealType<T>> void drawSphere(final RandomAccessibleInterval<T> img) {

		// the number of dimensions
		long[] center = new long[img.numDimensions()];
		for (int d = 0; d < img.numDimensions(); ++d) {
			center[d] = (img.max(d) - img.min(d)) / 2 + img.min(d);
		}

		// define the center and radius
		Point centerdef = new Point(center);
		long radius = 10;
		HyperSphere<T> testsphere = new HyperSphere<T>(img, centerdef, radius);

		Cursor<T> c = testsphere.cursor();

		long var = 0;

		while (c.hasNext()) {
			c.fwd();
			c.get().setReal(var);
			++var;

		}

	}

	public static void main(String[] args) {
		new ImageJ();

		final Img<FloatType> img = ArrayImgs.floats(1000, 1000);
		// drawSphere( img );
		// ImageJFunctions.show( img );
		long[] center = new long[img.numDimensions()];
		for (int d = 0; d < img.numDimensions(); ++d) {
			center[d] = (img.max(d) - img.min(d)) / 2 + img.min(d);
		}
		Point newcenter = new Point(center);
		drawMulti(img, newcenter, 400);
		ImageJFunctions.show(img);

	}

}
