Na pasta 003/:
	Modificar os caminhos no create_frames.py:
		IMAGE_PATH -> caminho para a pasta dos frames, dentro do 003/
		SNAPSGOT_PATH -> caminho para a pasta dos snapshots
	
	Modificar o caminho no join_frames.py:
		IMAGE_PATH -> caminho para a pasta dos frames, dentro do 003/
		PATH_TO_HDF5 -> caminho para a pasta dos snapshots
				
	Modificar o caminho no create_animation.sh:
		ls /mnt/d/IC/snapshots/snapshot_* > /mnt/d/IC/snapshots/snapshot.txt -> "/mnt/d/IC/snapshots/" deve ser o caminho até as snapshots

	Criar vídeo: ./create_animation.sh

Caso o create_animation.sh não esteja executável, escreva "chmod +x create_animation.sh" no terminal e execute
