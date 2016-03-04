package varun;



import net.imglib2.util.Pair;
import net.imglib2.util.ValuePair;

public class MergeSort2 {
	public static void sort() {

	}

	public static void sortMerge(int[] list, int[] listA, int[] listB) {

		int i=0, j=0, k=0;
		
		while(i<listA.length && j<listB.length){
			
			if(listA[i]<listB[j]){
				
				list[k]=listA[i];
				++i;
				++k;
			}
			
			
			else{
				
				list[k]=listB[j];
				++j;
				++k;
				
				
			}
			
			
		}
		
		while(i<listA.length){
			list[k]=listA[i];
			++i;
			++k;
			
		}
		
		while(j<listB.length){
			list[k]=listB[j];
			++j;
			++k;
			
		}
		
		
		
	}

	

	public static void split(int[] list) {

		if (list.length <= 1)
			return;
		else
		{

		int[] out1 = new int[list.length / 2];
		int[] out2 = new int[list.length / 2 + list.length % 2];

		int j = 0;
		for (int i = 0; i < out1.length; ++i) {
			out1[i] = list[j];
			++j;
		}

		for (int i = 0; i < out2.length; ++i) {
			out2[i] = list[j];
			++j;
		}

	//	Pair<int[], int[]> pair = new ValuePair<int[], int[]>(out1, out2);

		
			
			
			split(out1);
			
			
			split(out2);
			
			sortMerge(list,out1,out2);
			
			
			
		}
		
		}
		
	



	public static void main(String[] args) {

		int[] list = new int[8];
		

		list[0] = 5;
		list[1] = 2;
		list[2] = 4;
		list[3] = 7;
		list[4] = 1;
		list[5] = 3;
		list[6] = 2;
		list[7] = 6;

		split(list);
		
		for (int i=0; i<list.length; ++i){
			
			System.out.print(list[i]);
		}
		

	}

}
