public class Brain{

	private int neurons;
	private double[][] graph;
	
	public Brain(int neurons, double[][] graph){
		
		this.neurons = neurons;
		this.graph = graph;
		
	}
	
	public double[] brainMove(int[][] board){
		
		//everything except output
		double[] activity = new double[neurons - 4];
		
		for(int i = 0; i < 16; i++)
			activity[i] = board[i / 4][i % 4];
		
		for(int i = 16; i < neurons - 4; i++)
			for(int j = 0; j < neurons - 4; j++)
				activity[i] += activity[j]*graph[i][j];
		
		double[] output = new double[4];
		
		for(int i = 0; i < 4; i++)
			for(int j = 4; j < neurons - 4; j++)
				output[i] +=activity[j]*graph[i][j];
		
		return output;
		
	}
	
	public double[][] getGraph(){
		
		return graph;
		
	}
	
}
