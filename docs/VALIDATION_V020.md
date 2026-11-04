# Validação da infraestrutura de solvers v0.2.0 — 2026-11-04

## Preparação Git executada

Commit-base: **`e76c3197a73f43a00deeb7b9a0c4f08274df1e90`**,
`feat: validate geometry engine with real airfoils v0.1.1`.

A inspeção inicial encontrou `main` limpa, sem alterações ou arquivos relevantes
não rastreados. `origin` apontava para o repositório esperado:
`https://github.com/joaophuguet/AeroOpt.git`. Após `git fetch origin`, as referências
`main`, `origin/main`, `validation-v0.1.1` e `origin/validation-v0.1.1` continuavam
no mesmo commit. `git merge-base --is-ancestor validation-v0.1.1 main` retornou 0;
as comparações por `rev-list --left-right --count` retornaram `0 0`.

A v0.1.1 **já estava integralmente integrada e sincronizada**. Portanto não foi
necessário merge, pull nem push de `main`. A branch de destino não existia local
nem remotamente. Executou-se `git switch -c solver-interface-v0.2 main`; a nova
branch iniciou limpa e igual a `origin/main`, no commit-base acima.

Operações de preparação: `git status`, `git branch --show-current`, `git remote -v`,
`git branch -avv`, `git show`, `git fetch origin`, `git rev-list`,
`git merge-base --is-ancestor`, inspeção da ausência da nova branch e
`git switch -c solver-interface-v0.2 main`. Operações finais são inspeções de
status, diff, branch e commit. Não se alterou remoto, não se criou commit,
não se fez push, merge da v0.2.0, rebase, stash, limpeza ou exclusão de branch.
A autorização pontual de preparação não foi adicionada às regras permanentes.

## Ambiente e instalação

Windows 11, CPython 3.13.15. Todas as operações Python usaram
`C:\Users\stitc\Documents\Projetos\AeroOpt\.venv\Scripts\python.exe`.
Versão importada e versão instalada: **0.2.0**.

```powershell
$env:TEMP = Join-Path (Get-Location) '.tmp'
$env:TMP = $env:TEMP
$env:PIP_CACHE_DIR = Join-Path (Get-Location) '.cache\pip'
$env:PIP_DISABLE_PIP_VERSION_CHECK = '1'
.\.venv\Scripts\python.exe -m pip install -e '.[dev]' --no-deps
.\.venv\Scripts\python.exe -m pip check
```

Instalação editável concluída; `pip check`: `No broken requirements found.`
Nenhuma dependência de runtime nova, nenhuma instalação global, sem recriar
`.venv`. O controle Windows usa `ctypes`/Win32 e biblioteca padrão. NumPy 2.5.4,
Shapely 2.2.0 e pytest 9.2.0 continuam os mesmos. Caches, temporários, downloads
UIUC existentes e artefatos permanecem no workspace.

## Testes executados

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m pytest -q --reference-data=outputs/reference_data_v011
```

| Execução | Aprovados | Reprovados | Ignorados |
|---|---:|---:|---:|
| Padrão offline | 278 | 0 | 21 |
| Com os três originais UIUC preservados | **299** | **0** | **0** |

As últimas execuções completas terminaram em 25.35 s (offline) e 25.08 s (UIUC).

Os 21 ignorados no modo padrão são exclusivamente os testes publicados já
documentados na v0.1.1. Os **197 casos anteriores** continuam aprovados; foram
acrescentados **102 casos** de contratos, runner e processos. Não se alteraram
os módulos geométricos, suas tolerâncias, nem os arquivos de testes anteriores
para acomodar regressões. O modo com UIUC verifica os mesmos hashes fixados.

| Grupo | Verificações |
|---|---|
| Condições | Tipos, finitude, domínios básicos, ausência de limites arbitrários, Reynolds principal, alpha explícito, extensões tipadas |
| Imutabilidade | Condições, metadados aninhados, snapshots de coordenadas, preservação dos originais |
| Identidade | Ordenação canônica, mudanças de coordenadas/condições/solver/configuração/esquema/folga, independência de pasta e fonte, bits float64, zero com sinal |
| Resultados | Valores ausentes, rejeição de NaN/inf, referência de momento obrigatória, rótulo não físico, estados consistentes, reimportação e detecção de metadados alterados |
| Runner | Sete estados, geometria/folga inválida antes do preparo, falhas por etapa, repetição isolada, arquivos não sobrescritos, inventário e hashes |
| Processos | Saída zero/não zero, stdout/stderr, timeout real, encerramento de descendentes, falha de Job Object, saída normal com filho ainda ativo, ausência de shell |
| Auditoria | Saída malformada e arquivos parciais preservados, mutação de input rejeitada, SystemExit capturado, interrupção persistida e respeitada |
| Windows | Caminhos longos, espaços/Unicode, traversal, ADS, dispositivos, guarda de junction, Job Object real |
| Demonstração | Nove tentativas, todos os estados, mesma análise com UUID distinto, coeficientes ausentes, recusa de sobrescrita |

Arquivos novos de testes: `test_solver_contracts.py`, `test_solver_runner.py`,
`test_solver_processes.py`, com auxiliares em `solver_helpers.py` e
`fixtures/process_fixture.py`. Nenhum fixture implementa aerodinâmica. Um valor
artificial de coeficiente é usado apenas para verificar que sua publicação pelo
mock seja rejeitada; não aparece nos artefatos da demonstração.

## Evidência de timeout e encerramento

O fixture externo cria um programa-pai real e um filho real que escreve seu PID
e bloqueia por 3600 s. O pai confirma `CHILD_READY` somente após a criação do
filho. O runner aplica timeout de 1.5 s. Os testes verificam o estado `timeout`,
prazo global de teste inferior a 10 s, logs preservados e ambos os PIDs encerrados.
O Job Object inclui também o supervisor; seu contador precisa chegar a zero
antes de se declarar a limpeza completa.

Outro caso sai normalmente deixando o filho bloqueado. Mesmo com retorno zero,
o contexto encerra esse filho e verifica que ele não continua ativo. Também se
simula uma recusa de associação ao job: o alvo não chega a executar, não cria
arquivos de PID e o resultado é `execution_failed`, sem fallback sem contenção.

Há uma verificação separada em que o processo sai com código zero e o parser
declara convergência desconhecida: o estado é `not_converged`, não sucesso.
Um argumento contendo `&`, `>` e texto de comando é passado literalmente e
não cria o arquivo que uma interpretação por shell teria criado.

O teste de guarda de junction é feito por simulação do predicado do filesystem,
sem exigir privilégio de criação de symlink. Não se apresenta esse teste como
uma corrida hostil de filesystem. A contenção de processos foi exercitada
efetivamente no Windows; o backend POSIX existe mas não foi executado aqui.

## Problemas encontrados durante a implementação

Na primeira verificação do novo runner, arquivos temporários de publicação
JSON ultrapassaram o limite tradicional de caminho Windows, gerando
`FileNotFoundError` mesmo com a pasta existente. A correção usa o formato nativo
de caminhos longos, e há teste em uma pasta com caminho final acima de 260
caracteres, incluindo espaços e Unicode.

A primeira execução completa detectou três falhas nos testes de subprocesso:
o próprio Python externo não conseguia abrir `parent.pid` ou
`geometry/main.dat` por caminho relativo em um cwd muito longo. Suporte a
caminho longo no framework não garante suporte no programa externo. A estrutura
foi simplificada para `attempt-<UUID>` diretamente na raiz de execução, sem
concatenar o hash científico de 64 caracteres. O hash continua integralmente
preservado nos manifestos. Os 51 testes então existentes de runner/processos
passaram após essa correção, antes das verificações adicionais de tipos.

Outras proteções exercitadas: prazo/configuração inválidos não deixam NaN no
relatório; metadados são cópias profundas; strings não viram sequências de
caracteres por acidente; publicação de JSON não substitui destino existente;
`KeyboardInterrupt` salva evidências e é relançado; `SystemExit` de biblioteca
é convertido em falha estruturada. Nenhum teste antigo foi relaxado.

## Demonstração e inspeção dos artefatos

```powershell
.\.venv\Scripts\python.exe examples/demo_solver_interface.py --output outputs/solver_interface_v020_final
```

Entradas explícitas de demonstração: principal NACA 0012 analítico de 161
pontos; flap NACA 4412 analítico de 81 pontos; referências nominais por índices
explícitos; razão de corda 0.3; posição nominal `(1.02,-0.04)`; deflexão +25°;
alpha 6°, Re principal 200000 e Mach 0.05. Folga exigida 0.002 c; timeout 0.5 s.
São dados de teste da infraestrutura, sem recomendação de operação física.

| Caso | Resultado |
|---|---|
| success | success, convergência not_applicable |
| input_rejected | invalid_input |
| preparation_failed | preparation_failed, arquivo parcial preservado |
| execution_failed | execution_failed |
| timeout | timeout de processo externo real |
| parsing_failed | parsing_failed, saída malformada preservada |
| not_converged | not_converged |
| same_analysis_retry | success, mesmo case_id e UUID diferente |
| invalid_geometry | invalid_input antes de preparar o mock |

**Nove tentativas**: 2 sucessos simulados, 2 rejeições de entrada e um exemplar
de cada outro estado. Todos os resultados têm `cl=cd=cm=null`, identidade
simulada e `physical_coefficients_available=false`.

Foram reimportados os nove `result.json` e conferidos os hashes/tamanhos de
**81 artefatos**, a identidade derivada do snapshot e a separação dos UUIDs da
repetição. Logs do timeout registram `timed_out=true`, retorno 124 e
`windows_job_terminated_and_empty`. O prazo foi de 0.5 s, com pequeno tempo
adicional de limpeza e persistência (0.5105814 s registrados na evidência externa
final), sem confundir prazo de processo com
duração total da tentativa.

Saídas principais em `outputs/solver_interface_v020_final/`:

- `summary.json` com estados e caminhos relativos dos nove resultados.
- Uma pasta `attempt-<UUID>` por chamada, com `result.json`, eventos, inputs,
  DATs e arquivos das etapas alcançadas.
- `process.stdout.log`, `process.stderr.log`, `process.execution.json` e
  `execution.json` no cenário de timeout.
- `traceback.txt`, arquivo parcial ou saída inválida nas falhas correspondentes.

As execuções intermediárias em `outputs/solver_interface_v020/` e
`outputs/solver_interface_v020_verified/` foram preservadas. Para reproduzir,
escolha outra pasta relativa nova/vazia; não há sobrescrita silenciosa.
O teste da CLI também verifica explicitamente a recusa de uma segunda execução
na mesma pasta não vazia. Não foram gerados gráficos de sustentação/arrasto.

## Arquivos e decisões arquiteturais

Criados: pacote `src/aeroopt/solvers/` com os módulos públicos e auxiliares
documentados em `SOLVER_INTERFACE.md`; exemplo `demo_solver_interface.py`;
três arquivos de testes com seus dois auxiliares; `SOLVER_INTERFACE.md` e este
registro. Atualizados: `README.md`, `AGENTS.md`, `PROJECT_CONTEXT.md`, versão e
descrição em `pyproject.toml`, e `aeroopt.__version__`.

Decisões principais: dados imutáveis; snapshot global explícito; hash com tipo
e float64 hexadecimal; UUID separado; estados de processo e convergência
independentes; mock sem coeficientes; runner sequencial; logs em arquivo;
gate antes de associar/liberar processos; Job Object sem breakaway; caminhos
relativos contidos; JSON finito e publicação exclusiva por hard link em NTFS.
Não há dependência de solver, formato de malha ou ferramenta CFD específica.

## Limitações e prontidão

A infraestrutura passou pelos cenários exercitados e pode receber um adaptador
real em etapa posterior. Isso não valida física de solver nem classificação
de geometrias. O próximo adaptador precisa declarar seu binário/biblioteca e
versões, configurações cientificamente relevantes, faixas suportadas, condições
de transição, convenções de forças/momento e critérios de convergência lidos
da saída. Deve ser comparado com referências independentes e sensibilidade
à discretização antes de alimentar otimização.

Permanecem limites relevantes:

- Sem timeout preemptivo para código Python arbitrário em processo; bibliotecas
  bloqueantes devem usar o caminho externo para obter contenção verificável.
- O prazo não elimina latência de criação de processo no SO, preparação,
  parsing e I/O. Limpeza possui limite próprio.
- O contexto não é um sandbox para binários hostis ou escapes deliberados de
  job/grupo; adaptadores são confiáveis.
- Executáveis Windows legados podem precisar de raízes ainda mais curtas.
- Publicação JSON requer hard links; o conjunto de arquivos da tentativa não
  é transacional e falhas irrecuperáveis de disco/permissão podem impedir relatório.
- O hash não aproxima geometrias nem autentica autores; ambiente relevante e
  mudanças de código precisam ser declarados/versionados pelo adaptador.
- Não há promessa de resultados bit a bit entre versões de bibliotecas/SO.
- Gap/overlap aerodinâmicos, painelamento e toda a validação física continuam
  pendentes, conforme o contexto do projeto.
