# DeepSeek / OpenAI Function Calling Integration

## Option 1: Direct API

```python
import json
from openai import OpenAI  # or from deepseek import DeepSeek

client = OpenAI(api_key="...", base_url="https://api.deepseek.com/v1")

tools = json.load(open("../agent/tool_manifest.json"))["tools"]

response = client.chat.completions.create(
    model="deepseek-chat",
    messages=[
        {"role": "system", "content": open("../agent/system_prompt.md").read()},
        {"role": "user", "content": "Extract evidence spectra from these DIA files: ..."}
    ],
    tools=tools,
    tool_choice="auto"
)

# Handle tool calls
for choice in response.choices:
    if choice.finish_reason == "tool_calls":
        for call in choice.message.tool_calls:
            args = json.loads(call.function.arguments)
            if call.function.name == "run_extraction_chain":
                # Execute the chain
                result = run_chain(args)
                # Return result
                ...
```

## Option 2: LangChain

```python
from langchain.tools import tool
from langchain.agents import AgentExecutor, create_openai_tools_agent

@tool
def run_extraction_chain(search_results: str, raw_dir: str, output_dir: str,
                         protein_fasta: str = None, table1_xlsx: str = None) -> dict:
    """Run the full DIA spectrum evidence extraction chain."""
    import subprocess, sys
    # Execute scripts 70→71→73→76→78→79 in sequence
    scripts = [
        ("70_build_s14_targets.py", []),
        ("71_extract_targeted.py", []),
        ("73_profiles.py", []),
        ("76_build_isopsv_inputs.py", []),
        ("79_compose_a4.py", []),
    ]
    for script, extra in scripts:
        subprocess.run([sys.executable, f"scripts/{script}"] + extra, check=True)
    # R rendering
    subprocess.run(
        ["Rscript", "risoPSV/78_s14_export.R"],
        env={"OUT_PDF": "risoPSV/out/s14_topblock.pdf", "PATH": "..."},
        cwd="risoPSV"
    )
    return {"status": "done"}

agent = create_openai_tools_agent(llm, [run_extraction_chain], prompt)
executor = AgentExecutor(agent=agent, tools=[run_extraction_chain], verbose=True)
```

## Option 3: Direct CLI (no AI)

```bash
# Full chain
python scripts/70_build_s14_targets.py
python scripts/71_extract_targeted.py
python scripts/73_profiles.py
python scripts/76_build_isopsv_inputs.py
cd risoPSV && env -u LC_CTYPE Rscript 78_s14_export.R && cd ..
python scripts/79_compose_a4.py
python scripts/85_audit.py

# Table 1 cross-reference (requires extra chain)
python scripts/82_extra_targets.py
TGT=work/tmp/s14_targets_extra.parquet \
  OUTD=work/output/extract_extra \
  python scripts/71_extract_targeted.py
python scripts/84_extra_inputs.py
python scripts/81_table1_pdf.py
```

## Environment setup

| Variable | Description |
|----------|-------------|
| `DIA_EXTRACTION_ROOT` | Repository root (default: script parent dir) |
| `PEPTDEEP_CMD` | Path to peptdeep executable (default: `peptdeep`) |
| `RSCRIPT_CMD` | Path to Rscript (default: `Rscript`) |
