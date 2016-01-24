public class Main {

	public static void main(String[] args){
		
		Breeder breeder = new Breeder(10);

		while(breeder.maxFitness < 25000)
			breeder.breedNextGen();
		
	}
	
}
