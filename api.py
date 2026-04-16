import asyncio
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from combined import SUPPORTED_FORMATS, render_combined_figure

BASE_DIR = Path(os.getenv("RENDER_BASE_DIR", str(Path.cwd()))).resolve()
RENDER_TIMEOUT_SECONDS = int(os.getenv("RENDER_TIMEOUT_SECONDS", "120"))
MAX_CONCURRENT_RENDER = max(1, int(os.getenv("MAX_CONCURRENT_RENDER", "2")))
_render_semaphore = asyncio.Semaphore(MAX_CONCURRENT_RENDER)

app = FastAPI(title="50genes Visualization API", version="1.0.0")


class RenderRequest(BaseModel):
    csv_path: str = Field(..., description="CSV path (absolute or relative to base dir)")
    output_dir: str = Field(default="output", description="Output directory")
    filename_prefix: str = Field(default="Fig6", min_length=1)
    output_format: str = Field(default="svg", description="svg/png/pdf")


class RenderResponse(BaseModel):
    status: str
    output_file: str


def _resolve_safe_path(path_str: str) -> Path:
    path_obj = Path(path_str).expanduser()
    resolved = (BASE_DIR / path_obj).resolve() if not path_obj.is_absolute() else path_obj.resolve()
    if resolved != BASE_DIR and BASE_DIR not in resolved.parents:
        raise HTTPException(status_code=400, detail=f"Path escapes base directory: {path_str}")
    return resolved


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/render", response_model=RenderResponse)
async def render(req: RenderRequest) -> RenderResponse:
    output_format = req.output_format.lower()
    if output_format not in SUPPORTED_FORMATS:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported output_format: {output_format}. Supported: {sorted(SUPPORTED_FORMATS)}",
        )

    csv_path = _resolve_safe_path(req.csv_path)
    output_dir = _resolve_safe_path(req.output_dir)

    try:
        async with _render_semaphore:
            output_file = await asyncio.wait_for(
                asyncio.to_thread(
                    render_combined_figure,
                    csv_path=str(csv_path),
                    output_dir=str(output_dir),
                    filename_prefix=req.filename_prefix,
                    output_format=output_format,
                ),
                timeout=RENDER_TIMEOUT_SECONDS,
            )
        return RenderResponse(status="success", output_file=output_file)
    except TimeoutError:
        raise HTTPException(status_code=504, detail="Render timeout")
    except (FileNotFoundError, ValueError, PermissionError) as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except HTTPException:
        raise
    except Exception as exc:  # pragma: no cover
        raise HTTPException(status_code=500, detail=f"Unexpected server error: {exc}")


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("api:app", host="0.0.0.0", port=int(os.getenv("PORT", "8000")), reload=False)
