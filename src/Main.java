public class Main {

	public static void main(String[] args){
		
		Breeder breeder = new Breeder(100);
		System.out.println(breeder.maxFitness);
		while(breeder.generations < 100){
			breeder.breedNextGen();
			System.out.println(breeder.maxFitness);
		}
		
	}
	
}
