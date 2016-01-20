public class Main {

	public static void main(String[] args){
		
		Breeder breeder = new Breeder(100);

		while(breeder.generations < 100)
			breeder.breedNextGen();
		breeder.printFitness();
		
	}
	
}
