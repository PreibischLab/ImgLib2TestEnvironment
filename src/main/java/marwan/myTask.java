package marwan;

import marwan.Helper.Task;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class myTask extends Thread {

	private Portion currentPortion;
	private RandomAccessibleInterval<FloatType> targetView;
	private int num;
//	,sigma;
	private Task type ;
	
	

	public myTask(Portion portion, Img<FloatType> resultImage, Task type) {
		super();
		this.num = Helper.count++;
		this.currentPortion = portion;
		this.targetView= resultImage;
		this.type = type;
	}
	@Override
	public void run() {
		System.out.println("Thread"+num+" started.");
		switch (this.type) {
		case Gaus:
			try {
	Gauss3.gauss(Helper.sigma, this.currentPortion.getView(),Helper.targetPositon(targetView,currentPortion.getShape()));
	
			
			} catch (IncompatibleTypeException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;

		default:
			break;
		}
	
		System.out.println("Thread"+num+" finished.");
		if(Helper.log) ImageJFunctions.show(targetView).setTitle("target "+num);;
	}
}
