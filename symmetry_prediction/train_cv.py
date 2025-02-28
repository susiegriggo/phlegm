# imports 
import click 
from loguru import logger
import sys
from src import model
import os

logger.add("training.log")

# click stuff 
@click.command()
@click.option('--af2_embeddings_dir', type=click.Path(exists=True, file_okay=False), required=True, help='Directory containing AF2 embeddings')
@click.option('--output_dir', type=click.Path(file_okay=False), default='out', help='Directory to save output')
@click.option('--symmetries', type=click.Path(exists=True, file_okay=False), required=True, help='Path to symmetries file')
@click.option('--k_folds', type=int, default=5, help='Number of folds for cross-validation')
@click.option('--num_structures', type=int, default=5, help='Number of structures included in the AlphaFold embeddings')
@click.option('--use_pairwise', type=bool, default=True, help='Whether to include pairwise embeddings')
@click.option('--model_type', type=str, default='logistic_regression', help='Type of model to train')

def main(af2_embeddings_dir, output_dir, symmetries, k_folds, num_structures, use_pairwise, model_type): 
    try: 
        logger.info('Starting training')

        # Check if symmetries file is tab-separated
        with open(symmetries, 'r') as file:
            header = file.readline()
            if '\t' not in header:
                raise ValueError('Symmetries file must be tab-separated')

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Train the model
        model, results, internal_representations = model.train_model(
            af2_embeddings=af2_embeddings_dir,
            symmetries=symmetries,
            k_folds=k_folds,
            num_structures=num_structures,
            use_pairwise=use_pairwise,
            output_dir=output_dir,
            model_type=model_type
        )

        logger.info('Training completed successfully')

    except Exception as e: 
        logger.error(e)
        sys.exit(1) 
    
if __name__ == '__main__': 
    main()