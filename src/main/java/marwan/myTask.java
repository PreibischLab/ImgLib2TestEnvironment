package marwan;

import marwan.Helper.Task;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class myTask extends Thread {

	private Portion currentPortion;
	private RandomAccessibleInterval<FloatType> targetView;
	private int num;
	private Task type;

	public myTask(Portion portion, Img<FloatType> resultImage, Task type) {
		super();
		this.num = Helper.count++;
		this.currentPortion = portion;
		this.targetView = resultImage;
		this.type = type;
	}

	@Override
	public void run() {
		System.out.println("Thread" + num + " started.");
		switch (this.type) {
		case Gaus:
			try {
				if (this.currentPortion.getDimenssion() >0) {
					Gauss3.gauss(Helper.sigma, Views.extendMirrorSingle(this.currentPortion.getView()),
							Views.hyperSlice(targetView, this.currentPortion.getDimenssion(),
									this.currentPortion.getSlice()));
				} else {
					Gauss3.gauss(Helper.sigma, Views.extendMirrorSingle(this.currentPortion.getView()),
							targetView);
				}

			} catch (IncompatibleTypeException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;

		default:
			break;
		}

		System.out.println("Thread" + num + " finished.");
	}
}
