package varun;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;

import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;

import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class distanceTransform {

	public static void createBitimage(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<BitType> imgout,
			FloatType ThresholdValue) {

		final Cursor<FloatType> bound = Views.iterable(img).localizingCursor();

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

	public static <T extends RealType<T>> void computeDistance(RandomAccessibleInterval<BitType> img,
			RandomAccessibleInterval<T> imgout) {

		final Cursor<BitType> bound = Views.iterable(img).cursor();

		final RandomAccess<T> outbound = imgout.randomAccess();

		while (bound.hasNext()) {
			bound.fwd();
			outbound.setPosition(bound);

		
			if (bound.get().getInteger() == 0) {
				final Cursor<BitType> second = Views.iterable(img).cursor();
				double mindistance = 1.0E300;
				

				double distance = 0;
				while (second.hasNext())
					if (second.next().getInteger() == 1)

						for (int d = 0; d < second.numDimensions(); ++d) {

							distance += Math.pow((second.getDoublePosition(d) - bound.getDoublePosition(d)), 2);

						}

				mindistance = Math.min(mindistance, Math.sqrt(distance));
				
				System.out.println(mindistance);

				outbound.get().setReal(mindistance);
			}

			else {
				outbound.get().setReal(0);
			}

		}

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));

		ImageJFunctions.show(img).setTitle("Original_Image");

		final Img<BitType> imgout = new ArrayImgFactory<BitType>().create(img, new BitType());

		final Img<BitType> bitimgout = new ArrayImgFactory<BitType>().create(img, new BitType());

		FloatType val = new FloatType(200);

		createBitimage(img, imgout, val);

		computeDistance(imgout, img);

		ImageJFunctions.show(img).setTitle("FloatType_output");

		computeDistance(imgout, bitimgout);

		ImageJFunctions.show(bitimgout).setTitle("BitType_output");

	}

}
