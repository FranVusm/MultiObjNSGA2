# NSGA-II para TD-MT-GVRP biobjetivo

Este repositorio contiene una adaptacion de NSGA-II para resolver de forma heuristica un problema de ruteo de vehiculos verde, biobjetivo, dependiente del tiempo y con multiples viajes (TD-MT-GVRP).

El algoritmo minimiza simultaneamente:

- `F1`: emisiones de CO2.
- `F2`: costos operativos.

La implementacion parte de una plantilla clasica de NSGA-II, pero reemplaza la codificacion binaria original por una representacion especifica para rutas.

## Como compilar

En Linux, macOS o entornos con `make`:

```bash
make
```

Esto genera el ejecutable:

```bash
./nsga2r
```

En Windows, si ya existe el ejecutable compilado:

```powershell
.\nsga2r_port.exe
```

Si necesitas recompilar en Windows con GCC/MinGW, puedes usar el `Makefile` desde un entorno compatible con `make`, o compilar todos los `.c` enlazando la libreria matematica:

```bash
gcc *.c -o nsga2r_port.exe -lm
```

## Como ejecutar

La forma general de ejecucion es:

```bash
./nsga2r seed instance_route popsize ngen pcross pmut [runs]
```

En Windows:

```powershell
.\nsga2r_port.exe seed instance_route popsize ngen pcross pmut [runs]
```

Argumentos:

- `seed`: semilla base en el intervalo `(0,1)`.
- `instance_route`: ruta al archivo `.dat` de la instancia AMPL.
- `popsize`: tamano de poblacion. Debe ser mayor o igual a 4 y multiplo de 4.
- `ngen`: numero de generaciones.
- `pcross`: probabilidad de cruzamiento.
- `pmut`: probabilidad de mutacion.
- `runs`: numero opcional de corridas con semillas distintas. Si se omite, se usan 10.

Ejemplo:

```powershell
.\nsga2r_port.exe 0.123 ..\..\data\Modelo_juguete.dat 40 100 0.9 0.2 10
```

Otro ejemplo para una instancia mayor:

```powershell
.\nsga2r_port.exe 0.123 ..\..\data\Modelo_diverso.dat 32 600 0.8 0.3 10
```

## Datos de entrada

El lector espera un archivo `.dat` con formato AMPL. La instancia debe incluir, al menos:

- `N`: conjunto de nodos.
- `P`: conjunto de periodos.
- `K`: conjunto de vehiculos.
- `V`: conjunto de vueltas.
- `O`: deposito origen.
- `q`: capacidad del vehiculo.
- `d`: demanda por nodo.
- `LI`, `LS`: limites inferior y superior de cada periodo.
- `e`: emisiones por arco y periodo.
- `g`: costos por arco y periodo.
- `T`: tiempos de viaje por arco y periodo.
- `ee`: emisiones por detencion/servicio.
- `gg`: costos por detencion/servicio.
- `tt`: tiempos de servicio.

El ultimo nodo de `N` se interpreta como deposito ficticio de retorno.

## Salidas

Cada ejecucion crea una subcarpeta dentro de:

```text
resultados/
```

El nombre de la carpeta incluye el nombre de la instancia y un timestamp. Dentro se generan archivos como:

- `params.out`: parametros leidos y configuracion de ejecucion.
- `seed_XX_initial_pop.out`: poblacion inicial de cada semilla.
- `seed_XX_final_pop.out`: poblacion final de cada semilla.
- `seed_XX_best_pop.out`: soluciones factibles y no dominadas por semilla.
- `seed_XX_all_pop.out`: historial de poblaciones.
- `multi_seed_run_summary.csv`: resumen por semilla.
- `aggregate_pareto_solutions.txt`: soluciones agregadas no dominadas.
- `aggregate_pareto_points.csv`: puntos del frente agregado.
- `pareto_front.png`: grafico del frente agregado.

La carpeta `resultados/` se mantiene en Git con `.gitkeep`, pero sus datos de ejecucion estan ignorados para no subir resultados pesados o locales.

## Como funciona el algoritmo

NSGA-II mantiene una poblacion de soluciones candidatas. En cada generacion:

1. Evalua emisiones y costos de cada individuo.
2. Ordena las soluciones por dominancia de Pareto.
3. Usa distancia de crowding para preservar diversidad.
4. Selecciona padres por torneo binario.
5. Aplica cruzamiento y mutacion.
6. Mezcla padres e hijos y conserva las mejores soluciones no dominadas.

La ejecucion multi-semilla repite este proceso varias veces y luego agrega las soluciones de rango 1 para construir un frente Pareto global.

## Representacion del individuo

La adaptacion principal esta en la representacion. En vez de usar directamente variables binarias del tipo `x_ij_pvk`, cada individuo se codifica con:

```text
perm  = orden global de clientes
cuts  = separacion de ese orden en viajes/vehiculos
alpha = espera antes de iniciar cada viaje
```

### `perm`

`perm` es una permutacion de todos los clientes con demanda positiva. Cada cliente aparece exactamente una vez. Esto asegura que cada cliente sea atendido una sola vez.

Ejemplo:

```text
perm = [4, 2, 6, 3, 5]
```

### `cuts`

`cuts` divide la permutacion en segmentos. Cada segmento se asigna a un slot `(v,k)`, donde `v` es la vuelta y `k` el vehiculo.

Si hay `n_slots`, entonces:

```text
cantidad de cuts = n_slots + 1
n_slots = cantidad de vehiculos * cantidad de vueltas
```

Ejemplo:

```text
perm = [4, 2, 6, 3, 5]
cuts = [0, 2, 3, 5]
```

Esto genera:

```text
slot 0: [4, 2]
slot 1: [6]
slot 2: [3, 5]
```

Si dos cortes consecutivos son iguales, el slot queda vacio.

### `alpha`

`alpha` guarda la espera antes de iniciar cada slot. Permite que el algoritmo desplace una ruta hacia otro periodo temporal y aproveche diferencias de congestion, costos o emisiones.

## Decodificacion y evaluacion

Para evaluar un individuo, el codigo transforma cada segmento activo en una ruta:

```text
O -> clientes del segmento -> deposito ficticio
```

Luego recorre cada arco y calcula:

- periodo temporal de salida,
- emision del arco,
- costo del arco,
- tiempo de viaje,
- tiempo de servicio,
- hora acumulada para el siguiente arco.

Los objetivos se acumulan como:

```text
F1 += e[from,to,period] + ee[to,period]
F2 += g[from,to,period] + gg[to,period]
```

Si una solucion viola capacidad, estructura de rutas, horizonte temporal o periodos validos, recibe una penalizacion en `constr_violation`. Las soluciones factibles tienen `constr_violation = 0`.

## Operadores geneticos

El cruzamiento y la mutacion fueron adaptados a la representacion:

- Cruzamiento ordenado sobre `perm`, para evitar clientes repetidos o faltantes.
- Mezcla hereditaria de `cuts`, seguida de reparacion por capacidad.
- Promedio perturbado de `alpha`, acotado por limites temporales.
- Mutacion por intercambio de clientes.
- Mutacion por insercion de clientes.
- Mutacion por desplazamiento o reinicio de cortes.
- Mutacion de esperas temporales.

Despues de cruzar o mutar, el algoritmo repara los cortes cuando es necesario para mantener rutas coherentes y respetar capacidad.

## Scripts auxiliares

Para graficar un frente ya generado:

```powershell
python plot_pareto.py resultados\nombre_de_la_carpeta
```

Para probar combinaciones de parametros:

```powershell
python tune_parametros.py --model ..\..\data\Modelo_diverso.dat --runs 10 --max-attempts 20
```

El script de tuning registra combinaciones ya usadas en:

```text
tuning_combinaciones_usadas.txt
```

## Archivos principales

- `nsga2r.c`: rutina principal, ejecucion multi-semilla y reportes agregados.
- `reader.c`: lectura de instancias AMPL.
- `global.h`: estructuras globales, individuo y datos del problema.
- `initialize.c`: inicializacion y reparacion de cromosomas.
- `eval.c`: decodificacion, evaluacion y factibilidad.
- `crossover.c`: cruzamiento adaptado.
- `mutation.c`: mutaciones adaptadas.
- `dominance.c`: comparacion por factibilidad y dominancia Pareto.
- `fillnds.c`: ordenamiento no dominado y supervivencia.
- `report.c`: escritura de poblaciones.
- `plot_pareto.py`: grafico del frente agregado.
- `tune_parametros.py`: busqueda de parametros.
