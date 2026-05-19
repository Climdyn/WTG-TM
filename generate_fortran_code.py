
from model_definition import define_model

# constructing the model
model_definition = define_model(2, 2)

# computing the tensor (might take a long time depending on the resolution
model_definition.compute_tensor(numerical=False, compute_inner_products=True, compute_inner_products_kwargs={'timeout': None})

# generating the tendencies
f, Df = model_definition.compute_tendencies(language='fortran')

# writing to file
